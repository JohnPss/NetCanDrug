# Carregar pacote básico de sobrevivência (já vem no R)
library(survival)

# 1. Ler os dados (gerados no passo anterior)
data <- read.csv("survival_data_TOP2A.csv")

# 2. Preparar os dados
# Criar coluna de Tempo (escolhe death ou follow-up)
data$time <- ifelse(!is.na(data$days_to_death), data$days_to_death, data$days_to_last_follow_up)
# Criar coluna de Status (1=Morto, 0=Vivo)
# Verifica se o status contem "Dead" ou "Deceased"
data$status <- ifelse(grepl("Dead|Deceased", data$vital_status, ignore.case = TRUE), 1, 0)

# Filtrar dados válidos
data <- data[!is.na(data$time) & data$time > 0, ]

# Converter para Meses (Fica mais elegante no eixo X)
data$time_months <- data$time / 30.44

# Definir Grupos (Mediana)
cutoff <- median(data$expression, na.rm = TRUE)
data$group <- ifelse(data$expression > cutoff, "Alta Expressão (Tumor Agressivo)", "Baixa Expressão")
data$group <- factor(data$group, levels = c("Baixa Expressão", "Alta Expressão (Tumor Agressivo)"))

# 3. Calcular Estatísticas
surv_obj <- Surv(time = data$time_months, event = data$status)
fit <- survfit(surv_obj ~ group, data = data)

# Teste Log-Rank para pegar o P-valor
diff <- survdiff(surv_obj ~ group, data = data)
p_val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 4)))

# 4. Gerar Gráfico "Publication Ready" (Base R)
# Configuramos um PNG de alta resolução
png("results/figures/KM_TOP2A_NatureStyle.png", width = 2400, height = 2000, res = 300)

# Margens ajustadas
par(mar = c(6, 6, 4, 2) + 0.1)

# Plot principal
plot(fit,
    col = c("#2E9FDF", "#E7B800"), # Azul e Dourado (Padrão Nature/JCO)
    lwd = 4, # Linhas grossas (importante!)
    xlab = "Tempo Global de Sobrevida (Meses)",
    ylab = "Probabilidade de Sobrevida",
    main = "Impacto Prognóstico do gene TOP2A (TCGA-BRCA)",
    cex.lab = 1.5, # Tamanho labels eixos
    cex.axis = 1.3, # Tamanho numeros eixos
    cex.main = 1.6, # Tamanho titulo
    frame.plot = FALSE, # Remove a caixa em volta (estilo limpo)
    mark.time = TRUE # Marca os censurados com tracinhos
)

# Adicionar Eixos estilizados (Eixo L em vez de caixa)
axis(1, lwd = 2, cex.axis = 1.3)
axis(2, lwd = 2, cex.axis = 1.3)

# Adicionar Grid bem suave (cinza claro pontilhado)
grid(col = "gray90", lty = "dotted", lwd = 1)

# Adicionar Legenda
legend("bottomleft",
    legend = c("Baixa Expressão", "Alta Expressão (Pior Prognóstico)"),
    col = c("#2E9FDF", "#E7B800"),
    lwd = 4,
    bty = "n", # Sem caixa na legenda
    cex = 1.2
)

# Adicionar P-valor no gráfico
text(
    x = max(data$time_months) * 0.7, y = 0.85,
    labels = paste("Log-rank Test\n", p_text),
    cex = 1.4, font = 2, pos = 4
)

# Fechar e salvar
dev.off()

message("✅ Gráfico salvo em: results/figures/KM_TOP2A_NatureStyle.png")
