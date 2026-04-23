library(haven)
library(gtsummary)
library(pROC)
library(rsample)
library(ggplot2)
library(rms)
library(predtools)
library(readxl)
library(dcurves)
library(survival)

df <- read_xlsx("C:\\R\\Polyp\\data\\updated_polyp1.xlsx")

{
  df$Sex <- factor(df$Sex, levels = c(0,1), labels = c(0,1))
  
  df$Tuoicao <- ifelse(df$Age >= 50, 1, 0)
  df$Tuoicao <- factor(df$Tuoicao, levels = c(0,1), labels = c(0,1))
  
  df$HTL <- ifelse(df$Thuocla >= 1, 1, 0)
  df$HTL <- factor(df$HTL, levels = c(0,1), labels = c(0,1))
  
  df$rightcolon <- ifelse(df$Location >= 6, 1, 0)
  df$rightcolon <- factor(df$rightcolon, levels = c(0,1), labels = c(0,1))
  
  df$Paris <- factor(df$Paris, levels = c(0,1), labels = c(0,1))
  
  df$kichthuoclon <- ifelse(df$Size >= 6, 1, 0)
  df$kichthuoclon <- factor(df$kichthuoclon, levels = c(0,1), labels = c(0,1))
  
  df$WLE <- ifelse(df$WLE == 2, 1, 0)
  df$WLE <- factor(df$WLE, levels = c(0,1), labels = c(0,1))
  
  df$nhompolyp <- ifelse(df$nhompolyp == 1, 1, 0)
  df$nhompolyp <- factor(df$nhompolyp, levels = c(0,1), labels = c(0,1))
  
  df$TSUTDT <- factor(df$TSUTDT, levels = c(0,1), labels = c(0,1))
  
  df$TSGD <- factor(df$TSGD, levels = c(0,1), labels = c(0,1))
  
  df$Ruou <- factor(df$Ruou, levels = c(0,1), labels = c(0,1))
  
  df$BMI <- (df$CN / (df$Cao/100)^2)
  
  df$NBI <- ifelse(df$NBI == 2, 1, 0)
  df$NBI <- factor(df$NBI, levels = c(0,1), labels = c(0,1))
}


## Tổng quan dữ liệu hiện có
tbl_summary(df, by = nhompolyp, include = c(Age, Sex, TSUTDT, TSPDT, TSGD, HTL, Ruou, BMI, kichthuoclon, Paris, rightcolon, NBI),
            digits = c(all_categorical() ~ c(0, 1),
                       all_continuous() ~ 1),
            statistic = list(all_continuous() ~ "{mean} ({sd})")) %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3),
        test = list(all_continuous() ~ "t.test", all_categorical() ~ "chisq.test"))

## Phân chia dữ liệu
set.seed(123)
dfsplit <- initial_split(df, prop = 0.7)
train <- training(dfsplit)
test  <- testing(dfsplit)

## Vẽ đường cong ROC
# WLE
mod1 <- glm(nhompolyp ~ WLE, data = train, family = "binomial")
pred <- predict(mod1, test, type = "response")
roc1 <- roc(test$nhompolyp, pred, ci = T)

# Tuoicao + Paris + kichthuoclon + rightcolon
mod2 <- glm(nhompolyp ~ Tuoicao + Paris + kichthuoclon + rightcolon, data = train, family = "binomial")
pred <- predict(mod2, test, type = "response")
roc2 <- roc(test$nhompolyp, pred, ci = T)

# WLE + Tuoicao + Paris + kichthuoclon + rightcolon
mod3 <- glm(nhompolyp ~ WLE + Tuoicao + Paris + kichthuoclon + rightcolon, data = train, family = "binomial")
pred <- predict(mod3, test, type = "response")
roc3 <- roc(test$nhompolyp, pred, ci = T)

# NBI + WLE + Tuoicao + Paris + kichthuoclon + rightcolon
mod4 <- glm(nhompolyp ~ WLE + Tuoicao + Paris + kichthuoclon + rightcolon + NBI, data = train, family = "binomial")
pred <- predict(mod4, test, type = "response")
roc4 <- roc(test$nhompolyp, pred, ci = T)

# NBI
mod5 <- glm(nhompolyp ~ NBI, data = train, family = "binomial")
pred <- predict(mod5, test, type = "response")
roc5 <- roc(test$nhompolyp, pred, ci = T)

# ROC curves
roc.list <- list(roc1, roc2, roc3, roc4, roc5)
ggroc(roc.list, aes = "color", size = 1, legacy.axes = T) +
  annotate("text", x = 0.2, y = 0.3, size = 6, hjust = 0,
           label = paste0("NSAST + 4 yếu tố = ", sprintf("%.3f", roc3$auc))) +
  annotate("text", x = 0.2, y = 0.2, size = 6, hjust = 0,
           label = paste0("4 yếu tố = ", sprintf("%.3f", roc2$auc))) +
  annotate("text", x = 0.2, y = 0.1, size = 6, hjust = 0,
           label = paste0("NSAST = ", sprintf("%.3f", roc1$auc))) +
  annotate("text", x = 0.2, y = 0.4, size = 6, hjust = 0,
           label = paste0("NBI + 5 yếu tố = ", sprintf("%.3f", roc4$auc))) +
  annotate("text", x = 0.2, y = 0.5, size = 6, hjust = 0,
           label = paste0("NBI = ", sprintf("%.3f", roc5$auc))) +
  geom_point(x = 0.464, y = 0.82, color = "red", size = 4) +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  scale_color_discrete(name = "Mô hình", labels = c("NSAST", "4 yếu tố", "NSAST + 4 yếu tố", "NBI + 5 yếu tố", "NBI" )) +
  labs(x = "1 - Độ đặc hiệu", y = "Độ nhạy") + 
  ggtitle("Đường cong ROC của 5 mô hình")

## Độ nhạy - Độ đặc hiệu - Độ chính xác

# MODEL 1: WLE
coords(roc1, x = "best", best.method="youden", ret = c("threshold", "accuracy", "sensitivity", "specificity"))

# MODEL 2: 4 factors
coords(roc2, x = "best", best.method="youden", ret = c("threshold", "accuracy", "sensitivity", "specificity"))

# MODEL 3: WLE + 4 factors
coords(roc3, x = "best", best.method="youden", ret = c("threshold", "accuracy", "sensitivity", "specificity"))

# MODEL 4: NBI + 5 factors
coords(roc4, x = "best", best.method="youden", ret = c("threshold", "accuracy", "sensitivity", "specificity"))

# MODEL 5: NBI
coords(roc5, x = "best", best.method="youden", ret = c("threshold", "accuracy", "sensitivity", "specificity"))

## Mô hình NSAST + 4 yếu tố
auc(roc3)
ci.auc(roc3, conf.level=0.95, method=c("delong", "bootstrap"))
summary(mod3)


## Điểm cut-off của mô hình NSAST + 4 yếu tố
cutoff <- coords(roc3, x = "best", best.method="youden", ret = c("youden")) -1
print(paste0("Giá trị cut-off tối ưu trên ROC là ", round(cutoff,3)))

ggroc(roc3, aes = "color", size = 1, legacy.axes = T) +
  annotate("text", x = 0.4, y = 0.3, size = 6, hjust = 0,
           label = paste0("NSAST + 4 yếu tố = ", sprintf("%.3f", roc3$auc))) +
  geom_abline(color = "red") +
  theme_minimal() +
  theme(text = element_text(size = 18)) +
  labs(x = "1 - Độ đặc hiệu", y = "Độ nhạy") + 
  ggtitle("Đường cong ROC của mô hình NSAST + 4 yếu tố")


## Điểm C-index của mô hình
c_index <- concordance(mod3)
print(paste0("Điểm C-index là ", sprintf("%.3f", c_index["concordance"])))

## Điểm nguy cơ
round(coef(mod3) / min(abs(coef(mod3))), 0)

for (score in 0:17) {
  y <- as.numeric(min(abs(coef(mod3))) * score + coef(mod3)[1])
  p <- exp(y) / (1 + exp(y))
  print(paste0("Điểm ", score, " thì có nguy cơ p = ", sprintf("%.3f", p)))
}


## Biểu đồ nomogram
# Nomogram
d <- train[,c("nhompolyp", "WLE", "Tuoicao", "Paris", "kichthuoclon", "rightcolon")]
d <- zap_labels(d)
# Dán nhãn giá trị
nom_labels <- c(
  nhompolyp = "Loại Polyp",
  WLE = "Dự đoán trên NSAST",
  Tuoicao = "Tuổi",
  Paris = "Hình dạng đại thể",
  kichthuoclon = "Kích thước",
  rightcolon = "Vị trí"
)

label(d) <- lapply(names(nom_labels), function(x) label(d[,x]) <- nom_labels[x])
label(d)

ddist <- datadist(d)
options(datadist='ddist')

f <- lrm(nhompolyp ~ WLE + Tuoicao + Paris + kichthuoclon + rightcolon, data = d)

nom <- nomogram(f, fun = plogis, funlabel = "Risk", vnames = "labels")
plot(nom)

## Biểu đồ tỷ số số chênh
cols <- df |> select(Tuoicao, HTL, Ruou, TSGD, kichthuoclon, Paris, rightcolon, nhompolyp)

res.list <- lapply(as.list(cols), function(x) fisher.test(x, y = df$nhompolyp))
tmp <- do.call(rbind, lapply(res.list, broom::tidy))
tmp <- data.frame(tmp)

p <- lapply(as.list(cols), function(x) chisq.test(x, y = df$nhompolyp)$p.value)
p <- c(do.call(rbind, p))

tmp$feat <- colnames(cols)
tmp$feat <- c("Tuổi", "Hút thuốc lá", "Sử dụng rượu bia", "Tiền sử ung thư gia đình", "Kích thước", "Hình dạng đại thể", "Vị trí", Inf)
tmp$feat <- factor(tmp$feat, levels = rev(c("Tuổi", "Hút thuốc lá", "Sử dụng rượu bia", "Tiền sử ung thư gia đình", "Kích thước", "Hình dạng đại thể", "Vị trí")))
tmp$p.chisq <- p
tmp$p.chisq <- sprintf("%.3f", tmp$p.chisq)
tmp$p.chisq <- as.character(tmp$p.chisq)
tmp$p.chisq <- ifelse(tmp$p.chisq == "0.000", "<0.001", tmp$p.chisq)
tmp$range <- paste0(sprintf("%.1f", tmp$conf.low), "-", sprintf("%.1f", tmp$conf.high))

ggplot(tmp[1:7,], aes(x = feat, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 1.2) +
  geom_hline(yintercept = 1.0, linetype = "dashed", size = 1) +
  scale_y_log10(breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10), minor_breaks = NULL) +
  labs(y = "Odds ratio", x = NULL) +
  coord_flip(ylim = c(0.1, 20)) +
  theme_bw() +
  annotate("text", label = tmp$p.chisq, x = tmp$feat, y = 0.11) +
  annotate("text", label = tmp$range, x = tmp$feat, y = 15) +
  annotate("text", label = sprintf("%.2f", tmp$estimate), x = tmp$feat, y = tmp$estimate, vjust = -1) +
  theme(text = element_text(size = 15))

## Biểu đồ hiệu chuẩn
# Calibration plot
d <- train[,c("nhompolyp", "WLE", "Tuoicao", "Paris", "kichthuoclon", "rightcolon")]
f <- lrm(nhompolyp ~ WLE + Tuoicao + Paris + rightcolon + kichthuoclon, data = d)

d$prob_pred <- predict(f, type = "fitted")

d$prob_bins <- cut(d$prob_pred, breaks = seq(0, 1, by = 0.1), labels = FALSE)

d_calib <- aggregate(cbind(prob_pred, nhompolyp) ~ prob_bins, data = d, mean)

# Trừ 1 ở y do nhompolyp được mặc định là 1 và 2
ggplot(d_calib, aes(x = prob_pred, y = nhompolyp - 1)) +
  geom_point(size = 2.5, color = "darkgreen") +
  geom_line(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  # scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(x = "Mean Predicted Probability", y = "Observed Probability") +
  ggtitle("Biểu đồ hiệu chuẩn") +
  theme_bw()


## Phân tích đường cong quyết định DCA
# DCA nhóm test
test$validation <- predict(mod3, test, type = "response")
dcurve <- dca(nhompolyp ~ validation, 
              data = test,
              as_probability = c("validation"))

# DCA nhóm train
train$training <- predict(mod3, train, type = "response")
dcurve1 <- dca(nhompolyp ~ training,
              data = train,
              as_probability = c("training"))

plot(dcurve)
plot(dcurve1)

