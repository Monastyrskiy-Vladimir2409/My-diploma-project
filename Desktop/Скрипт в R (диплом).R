# -------------------- 1. Установка и загрузка пакетов --------------------
required_packages <- c("HardyWeinberg", "ggplot2", "logistf")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    library(pkg, character.only = TRUE)
  }
}

# -------------------- 2. Ввод данных --------------------
# Формат: c(гомозигота по референсному аллелю, гетерозигота, гомозигота по минорному аллелю)

PE_geno <- list(
  rs2854117 = c(13, 22, 7), # CC, CT, TT
  rs2854116 = c(11, 25, 6), # СС, ТС, ТТ (min=T)
  rs5128 = c(25, 14, 3), # CC, CG, GG
  rs5742904 = c(42, 0, 0), # GG, GA, AA
  rs17231506 = c(18, 17, 7), # CC, CT, TT
  rs328 = c(31, 11, 0), # CC, CG, GG
  rs9939609 = c(18, 14, 10) # TT, TA, AA
)

Control_geno <- list(
  rs2854117 = c(15, 27, 8),
  rs2854116 = c(12, 26, 12), # СС, ТС, ТТ
  rs5128 = c(43, 7, 0),
  rs5742904 = c(49, 1, 0),
  rs17231506 = c(26, 19, 5),
  rs328 = c(41, 9, 0),
  rs9939609 = c(19, 21, 10)
)

snp_names <- names(PE_geno)

# Минорные аллели (для справки, для расчета MAF)
minor_allele <- c(rs2854117="T", rs2854116="T", rs5128="G", rs5742904="A", rs17231506="T", rs328="G", rs9939609="A")

# -------------------- 3. Проверка равновесия Харди‑Вайнберга (HWE) --------------------
get_hwe_p <- function(obs) {
  if (all(obs[2:3] == 0)) return(NA)
  HWExact(obs, pvaluetype = "selome")$pval
}

hwe_table <- data.frame(
  SNP = rep(snp_names, 2),
  Group = rep(c("PE", "Control"), each = length(snp_names)),
  p_value = c(sapply(PE_geno, get_hwe_p), sapply(Control_geno, get_hwe_p))
)
print("HWE p‑значения:")
print(hwe_table)

# -------------------- 4. Ассоциация (АДДИТИВНАЯ модель) --------------------
# Строим датафрейм с индивидуальными данными для логистической регрессии
# Для этого преобразуем агрегированные генотипы в individual-level данные

# Функция для создания датафрейма из агрегированных генотипов
expand_geno <- function(geno_counts, group_label) {
  # geno_counts: c(гом_ref, гетеро, гом_мин)
  # Возвращает data.frame с колонками: Status (0/1), Add_score (0,1,2)
  n_ref <- geno_counts[1]
  n_het <- geno_counts[2]
  n_min <- geno_counts[3]
  status <- ifelse(group_label == "PE", 1, 0)
  data.frame(
    Status = rep(status, n_ref + n_het + n_min),
    Add_score = c(rep(0, n_ref), rep(1, n_het), rep(2, n_min))
  )
}

# Создаём единый датафрейм для всех SNP и групп
all_data_list <- list()
for (snp in snp_names) {
  pe_df <- expand_geno(PE_geno[[snp]], "PE")
  ctrl_df <- expand_geno(Control_geno[[snp]], "Control")
  df <- rbind(pe_df, ctrl_df)
  df$SNP <- snp
  all_data_list[[snp]] <- df
}
all_data <- do.call(rbind, all_data_list)

# Теперь для каждого SNP проводим логистическую регрессию Status ~ Add_score
assoc_additive <- data.frame()
for (snp in snp_names) {
  sub <- subset(all_data, SNP == snp)
  # Если в группе ПЭ или контроле все наблюдения имеют одинаковый Add_score (например, rs5742904 в ПЭ все 0)
  # то обычная логистическая регрессия может не сходиться. Используем поправку Фирта (logistf)
  if (var(sub$Add_score) == 0) {
    # все одинаковые – нет изменчивости, OR не определён
    OR <- NA
    CI_low <- NA
    CI_high <- NA
    p_val <- NA
  } else {
    # Проверяем на идеальное разделение (например, все носители минорного аллеля только в одной группе)
    # Используем logistf для устойчивости
    fit <- logistf(Status ~ Add_score, data = sub, pl = TRUE)
    OR <- exp(coef(fit)[2])
    CI_low <- exp(confint(fit)[2, 1])
    CI_high <- exp(confint(fit)[2, 2])
    p_val <- fit$prob[2]
  }
  assoc_additive <- rbind(assoc_additive,
                          data.frame(SNP = snp, OR = OR, CI_low = CI_low, CI_high = CI_high, P_value = p_val))
}

# Применяем поправку на множественные сравнения (Бенджамини-Хохберг)
assoc_additive$P_adj <- p.adjust(assoc_additive$P_value, method = "BH")
assoc_additive$OR_CI <- sprintf("%.2f (%.2f-%.2f)", assoc_additive$OR, assoc_additive$CI_low, assoc_additive$CI_high)

print("Ассоциация (аддитивная модель – OR на одну копию минорного аллеля):")
print(assoc_additive[, c("SNP", "OR", "CI_low", "CI_high", "P_value", "P_adj")])

# -------------------- 5. Сравнение с популяционными данными gnomAD --------------------
gnomad_maf <- c(
  rs2854117 = 0.29,
  rs2854116 = 0.375,
  rs5128 = 0.0876,
  rs5742904 = 0.00049,
  rs17231506= 0.31,
  rs328 = 0.089,
  rs9939609 = 0.39
)

compare_gnomad <- function(geno, n_people, maf_g) {
  minor_obs <- geno[2] + 2*geno[3]
  major_obs <- 2*n_people - minor_obs
  minor_exp <- 2*n_people * maf_g
  major_exp <- 2*n_people - min_exp
  tab <- matrix(c(minor_obs, major_obs, minor_exp, major_exp), nrow = 2, byrow = TRUE)
  fisher.test(tab)$p.value
}

gnomad_p_PE <- sapply(snp_names, function(s) compare_gnomad(PE_geno[[s]], 42, gnomad_maf[s]))
gnomad_p_Ctrl<- sapply(snp_names, function(s) compare_gnomad(Control_geno[[s]], 50, gnomad_maf[s]))

gnomad_res <- data.frame(
  SNP = snp_names,
  MAF_gnomAD = gnomad_maf,
  MAF_PE = sapply(snp_names, function(s) (PE_geno[[s]][2] + 2*PE_geno[[s]][3]) / 84),
  MAF_Ctrl = sapply(snp_names, function(s) (Control_geno[[s]][2] + 2*Control_geno[[s]][3]) / 100),
  p_vs_gnomad_PE = gnomad_p_PE,
  p_vs_gnomad_Ctrl = gnomad_p_Ctrl
)
print("Сравнение с gnomAD:")
print(gnomad_res)

# -------------------- 6. ГРАФИКИ --------------------

# 6.1 График forest.plot для OR (аддитивная модель)
assoc_plot <- assoc_additive[!is.na(assoc_additive$OR), ]
assoc_plot <- assoc_plot[order(assoc_plot$OR), ]
assoc_plot$SNP <- factor(assoc_plot$SNP, levels = assoc_plot$SNP)

forest_plot <- ggplot(assoc_plot, aes(x = OR, y = SNP)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, color = "blue") +
  geom_point(size = 3, color = "red") +
  scale_x_log10() +
  labs(x = "Odds ratio на одну копию минорного аллеля (95% CI)", y = "",
       title = "Ассоциация SNP с преэклампсией (аддитивная модель)") +
  theme_minimal()
print(forest_plot)
ggsave("Forest_plot_additive.png", width = 8, height = 5)

# 6.2 Столбчатая диаграмма частоты минорного аллеля (MAF) по группам
maf_df <- data.frame(
  SNP = rep(snp_names, 2),
  Group = rep(c("PE", "Control"), each = length(snp_names)),
  MAF = c(
    sapply(snp_names, function(s) (PE_geno[[s]][2] + 2*PE_geno[[s]][3]) / 84),
    sapply(snp_names, function(s) (Control_geno[[s]][2] + 2*Control_geno[[s]][3]) / 100)
  )
)

barplot_maf <- ggplot(maf_df, aes(x = SNP, y = MAF, fill = Group)) +
  geom_col(position = position_dodge(width = 0.9)) +
  labs(y = "Частота минорного аллеля (MAF)", title = "Сравнение MAF в группах") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(barplot_maf)
ggsave("MAF_comparison.png", width = 8, height = 5)

# 6.3 График «доза–ответ»: предсказанная вероятность ПЭ в зависимости от числа минорных аллелей
# Строим для SNP с наименьшим p‑значением (скорее всего rs5128)
best_snp <- assoc_additive$SNP[which.min(assoc_additive$P_value)]
if (!is.na(best_snp)) {
  sub_best <- subset(all_data, SNP == best_snp)
  # Повторная регрессия для этого SNP для получения предсказанных вероятностей
  fit_best <- glm(Status ~ Add_score, data = sub_best, family = binomial)
  newdata <- data.frame(Add_score = 0:2)
  pred <- predict(fit_best, newdata, se.fit = TRUE, type = "response")
  dose_df <- data.frame(Add_score = 0:2,
                        Prob = pred$fit,
                        Lower = pred$fit - 1.96 * pred$se.fit,
                        Upper = pred$fit + 1.96 * pred$se.fit)
  
  dose_plot <- ggplot(dose_df, aes(x = Add_score, y = Prob)) +
    geom_line(size = 1.2, color = "darkgreen") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2) +
    geom_point(size = 3) +
    scale_x_continuous(breaks = 0:2, labels = c("0", "1", "2")) +
    labs(x = "Число минорных аллелей", y = "Предсказанная вероятность ПЭ",
         title = paste("Доза-эффект для", best_snp, "(аддитивная модель)")) +
    theme_minimal()
  print(dose_plot)
  ggsave(paste0("Dose_response_", best_snp, ".png"), width = 6, height = 5)
} else {
  message("Не удалось определить SNP с минимальным p‑значением.")
}

# 6.4 График HWE: наблюдаемые vs ожидаемые генотипы (для всех SNP и групп)
expected_geno <- function(obs) {
  n <- sum(obs)
  p <- (2*obs[1] + obs[2]) / (2*n)
  q <- 1 - p
  expected <- c(p^2, 2*p*q, q^2) * n
  names(expected) <- names(obs)
  return(expected)
}

build_hwe_long <- function(geno_list, group_name) {
  do.call(rbind, lapply(snp_names, function(snp) {
    obs <- geno_list[[snp]]
    exp <- expected_geno(obs)
    geno_types <- c("Гом. реф.", "Гетеро.", "Гом. мин.")
    data.frame(SNP = snp, Group = group_name,
               Genotype = geno_types,
               Observed = obs, Expected = exp)
  }))
}

hwe_PE <- build_hwe_long(PE_geno, "ПЭ")
hwe_Ctrl <- build_hwe_long(Control_geno, "Контроль")
hwe_all <- rbind(hwe_PE, hwe_Ctrl)
hwe_long <- tidyr::pivot_longer(hwe_all, cols = c(Observed, Expected), names_to = "Type", values_to = "Count")

hwe_barplot <- ggplot(hwe_long, aes(x = Genotype, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  facet_grid(Group ~ SNP, scales = "free_y") +
  labs(title = "Равновесие Харди–Вайнберга: наблюдаемые и ожидаемые генотипы",
       y = "Число генотипов", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10, face = "bold"))
print(hwe_barplot)
ggsave("HWE_all_SNPs_bars.png", width = 14, height = 8)

# 6.5 Точечный график p‑значений HWE (с порогом 0,05)
hwe_pvals <- hwe_table
hwe_pvals$logp <- -log10(hwe_pvals$p_value)
hwe_pvals$Significant <- ifelse(!is.na(hwe_pvals$p_value) & hwe_pvals$p_value < 0.05, "p<0.05", "NS")

pval_plot <- ggplot(hwe_pvals, aes(x = SNP, y = logp, color = Group, shape = Significant)) +
  geom_point(size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(y = "-log10(p-value HWE)", title = "Проверка равновесия Харди–Вайнберга",
       subtitle = "Красная линия соответствует p=0,05") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("PE" = "blue", "Control" = "darkgreen")) +
  scale_shape_manual(values = c("p<0.05" = 17, "NS" = 16))
print(pval_plot)
ggsave("HWE_pvalues_points.png", width = 8, height = 5)
