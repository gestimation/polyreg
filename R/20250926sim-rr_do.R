library(mets)
library(prodlim)
data(bmt)



#cif1 <- cbind(c(0,10,20,100),c(0,0.15,0.2,0.22))
#fcif(cif1,c(1,0.5))

#cif1 <- cbind(c(0,10,20,100),c(0,0.15,0.2,0.22))
#cif2 <- cbind(c(0,10,20,100),c(0,0.1,0.15,0.18))
#cif2 <- cbind(c(0,10,20,100),c(0,0.3,0.35,0.38))
#cif2 <- cbind(c(0,10,20,100),c(0,0.5,0.55,0.58))
cif1 <- cbind(c(0,10,20,100),c(0,0.15,0.2,0.22))
cif2 <- cbind(c(0,10,20,100),c(0,0.7,0.75,0.78))

#simRR <- function(n,lrr1=c(-0.69314718056,-0.69314718056),lrr2=c(-0.69314718056,-0.69314718056),cens=NULL,
#simRR <- function(n,lrr1=c(0,-0.693147180560),lrr2=c(-0.69314718056,-0.69314718056),cens=NULL,
#simRR <- function(n,lrr1=c(-0.693147180560,-0.69314718056),lrr2=c(0,-0.69314718056),cens=NULL,
#simRR <- function(n,lrr1=c(0,-0.69314718056),lrr2=c(0,-0.69314718056),cens=NULL,
lrr1_10051005=c(0,-0.69314718056)
lrr2_10051005=c(0,-0.69314718056)
lrr1_05050505=c(-0.69314718056,-0.69314718056)
lrr2_05050505=c(-0.69314718056,-0.69314718056)
lrr1_10050505=c(0,-0.69314718056)
lrr2_10050505=c(-0.69314718056,-0.69314718056)
lrr1_05051005=c(-0.69314718056,-0.69314718056)
lrr2_05051005=c(0,-0.69314718056)

#rep=1000
n=20000
rep=1
simulation_size=rep*n

#type1_="cif"
#type2_="cif"
#type1_="logistic"
#type2_="cif"
#type1_="cloglog"
#type2_="cif"
type1_="cif"
type2_="cif"
cens_=0.01

T1 <- simRR(2000,lrr1_05050505,lrr2_05050505,cens_,type1_,type2_)
output <- polyreg(nuisance.model=Event(time, epsilon)~L, exposure = "A", data=T1, effect.measure1="RR", effect.measure2="RR", time.point=1, outcome.type = "COMPETING-RISK")
library(modelsummary)
msummary(output$summary, exponential=TRUE)

# --- 設定（必要に応じて編集） ---------------------------------------------
PATH1          <- "C:/path/to/dir"     # フォルダ
CSV_BASENAME   <- "input_name"         # 拡張子なしの CSV 名
out_csv        <- "res.csv"            # 出力CSVファイル名
rep            <- 100                  # 反復回数（&rep.）
size           <- 1000                 # 1反復あたりの行数（&size.）
EFFECT_MEASURE <- "RR"                 # &EFFECT_MEASURE.
TIMEPOINT      <- c(12, 24)            # &TIMEPOINT. 例：評価時点（月など）
S_NUM          <- 2                    # S_NUM=2（Sの水準数）
L_BIN          <- 1                    # L_BIN=1（Lがバイナリ等の指定）
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  if (requireNamespace("polyreg", quietly = TRUE)) library(polyreg)
  if (requireNamespace("dplyr", quietly = TRUE))   library(dplyr)
})
output <- polyreg(nuisance.model=Event(time, epsilon)~L, exposure = "A", data=T1, effect.measure1="RR", effect.measure2="RR", time.point=1, outcome.type = "COMPETING-RISK")
res <- output$summary$event1$tidy

# 解析1回分のフック関数：
# ここだけ、あなたの polyreg 呼び出し仕様に合わせて中身を書き換えてください。
run_polyreg <- function(data, effect_measure, timepoint, s_num = 2, l_bin = 1) {
  # data には列: id, time, epsilon, a, l（SASコード準拠。下で s = l + 1 を作ります）
  # 例：あなたの polyreg() が「Surv(time, epsilon) ~ a + s + L...」の形を取る場合は、
  # ↓の疑似コードを、実際の引数名・返り値に合わせて修正してください。

  # --- ダミー実装（返すのは data.frame）。必ず置き換えてください。 ---
  # 期待する返り値の形：少なくとも i_rep ごとの結果が1つの data.frame で返る。
  # 下は“仮”の出力（係数やCIがあればそれを入れる）。失敗時は NA を返す。
  out <- data.frame(
    id_min        = min(data$id),
    id_max        = max(data$id),
    n_rows        = nrow(data),
    effect_measure= effect_measure,
    timepoint     = paste(timepoint, collapse = ";"),
    estimate      = NA_real_,
    conf_low      = NA_real_,
    conf_high     = NA_real_
  )
  return(out)
  # -------------------------------------------------------------------------

  # 参考：もし polyreg::polyreg() が次のように呼べるなら（仮例）
  # fit <- polyreg::polyreg(
  #   formula = survival::Surv(time, epsilon) ~ a + s,   # L は L_LIST で別指定等なら調整
  #   data    = data,
  #   effect.measure = effect_measure,
  #   time.point     = timepoint,
  #   s.num          = s_num,
  #   l.bin          = l_bin,
  #   id.var         = "id"            # ID を渡す必要があれば
  # )
  # # 返却：係数表などを data.frame 化
  # broom::tidy(fit) |>
  #   dplyr::mutate(effect_measure = effect_measure,
  #                 timepoint = paste(timepoint, collapse = ";"))
}

# --- 入力を読み込み（SAS: PROC IMPORT の代替） ----------------------------
infile <- file.path(PATH1, paste0(CSV_BASENAME, ".csv"))
if (!file.exists(infile)) stop("入力CSVが見つかりません: ", infile)

data_raw <- utils::read.csv(infile, check.names = FALSE)
# 行ID付与（SAS: id=_N_;）
data_raw$id <- seq_len(nrow(data_raw))


size <- 2000
rep <- 1
i_rep <- 1
data_raw <- simRR(size*rep,lrr1_05050505,lrr2_05050505,cens_,type1_,type2_)
data_raw$id <- seq_len(nrow(data_raw))
res <- NULL
for (i_rep in seq_len(rep)) {
  id_min <- (i_rep - 1L) * size + 1L
  id_max <- i_rep * size
  chunk  <- data_raw[data_raw$id >= id_min & data_raw$id <= id_max, , drop = FALSE]
  output <- try(
    polyreg(nuisance.model=Event(time, epsilon)~L, exposure="A", data=chunk, effect.measure1="RR", effect.measure2="RR", time.point=1, outcome.type = "COMPETING-RISK"),
    silent = TRUE
  )
  if (inherits(output, "try-error") || is.null(output)) {
    warn_msg <- if (inherits(output, "try-error")) as.character(output) else "NULL returned"
    res_tmp <- data.frame(
      id_min         = id_min,
      id_max         = id_max - 1L,
      n_rows         = nrow(chunk),
      effect_measure = "RR",
      timepoint      = 1,
      estimate       = output$summary$event1$tidy$estimate,
      conf_low       = output$summary$event1$tidy$conf_low,
      conf_high      = output$summary$event1$tidy$conf_high,
      error          = paste0("i_rep=", i_rep, ": ", warn_msg),
      stringsAsFactors = FALSE
    )
  } else {
    res_tmp$i_rep <- i_rep
  }
  res <- if (is.null(res)) res_tmp else dplyr::bind_rows(res, res_tmp)
}




utils::write.csv(res, file = file.path(PATH1, out_csv), row.names = FALSE)
message("書き出し完了: ", file.path(PATH1, out_csv))
