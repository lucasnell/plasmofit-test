## How long is the ring stage in 3D7, and what does that make the inoculum's
## age window? (2026-11-09, Lucas, after he supplied llinas-2006.pdf.)
##
## TWO INDEPENDENT SOURCES, both microscopy staging of synchronized cultures:
##  Smith et al. (2020) Science, Table S6: a FITTED ring-to-trophozoite
##    transition phase per strain (theta_R,T). 3D7 = 0.312 of the cycle.
##    Their Table S5 gives the same as 12 h of a 39 h period (30.8%).
##  Llinas et al. (2006) NAR 34:1166, Figure 1: hourly staging over the IDC
##    for 3D7, Dd2 and HB3, drawn as a stacked grey band above each
##    phaseogram. There is no table, so the percentages are READ OFF THE
##    FIGURE by pixel classification -- `llinas-2006.pdf` rendered at 400 dpi,
##    the 3D7 band cropped, each column classified by grey level into ring
##    (darkest), trophozoite, schizont, with the intermediate grey at the
##    ring/trophozoite boundary counted as trophozoite. This is a measurement
##    off a figure and should be quoted as approximate.
##
## The time axis: Llinas measured 53 h for 3D7, and the next cycle's rings
## reappear partway along the panel, so the IDC is taken as the span from the
## panel's left edge to the reappearance of rings, and boundaries are reported
## as a FRACTION of that, then scaled to a 48 h cycle.
##
## Why it matters: the IBSM inoculum is ring-selected (circulating-only plus
## Glycerolyte freezing), so the ring window, not the wider circulating
## window, is the support the starting age distribution should have
## (`claude/start-dist-design.md`).
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/ring-stage-window.R > _data/ring-stage-window-<date>.txt
##   (needs llinas-2006.pdf in the repo root; pdftoppm and ImageMagick convert)

PDF <- "llinas-2006.pdf"
TMP <- file.path(Sys.getenv("CLAUDE_JOB_DIR", tempdir()), "tmp", "ringstage")
dir.create(TMP, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Smith et al., the fitted transition phases ------------------------
SMITH <- data.frame(
    strain = c("3D7", "HB3", "SA250", "D6", "FVO-NIH"),
    period = c(39, 46, 54, 36, 43),
    theta_RT = c(0.312, 0.390, 0.419, 0.478, 0.571))
SMITH$ring_h_own <- round(SMITH$theta_RT * SMITH$period, 1)
SMITH$ring_h_48 <- round(SMITH$theta_RT * 48, 1)
cat("1. SMITH ET AL. (2020) Table S6: fitted ring-to-trophozoite transition\n",
    "   phase, as a fraction of the cycle, at each strain's own period and\n",
    "   rescaled to a 48 h cycle.\n\n", sep = "")
print(SMITH, row.names = FALSE)

## ---- 2. Llinas et al., read off Figure 1 ----------------------------------
png <- file.path(TMP, "p4")
system2("pdftoppm", c("-f", "4", "-l", "4", "-r", "400", "-png", PDF, png))
p4 <- paste0(png, "-4.png")
stopifnot(file.exists(p4))
txt <- file.path(TMP, "px.txt")
system(sprintf("convert %s -crop 340x135+548+500 +repage -colorspace gray txt:- | tail -n +2 | sed 's/[,:]/ /g' | awk '{print $1, $2, $3}' > %s",
               shQuote(p4), shQuote(txt)))
d <- read.table(txt, sep = " ", stringsAsFactors = FALSE)
d$x <- as.integer(d[[1]]); d$y <- as.integer(d[[2]])
d$g <- as.numeric(gsub("[()]", "", d[[3]]))
d <- subset(d, y >= 10 & y <= 132)                 # the plot area, 123 px = 100%
b <- ifelse(d$g < 95, "ring", ifelse(d$g < 150, "mid",
     ifelse(d$g < 190, "troph", ifelse(d$g < 240, "schiz", "white"))))
t <- as.data.frame.matrix(table(d$x, b)); t$x <- as.integer(rownames(t))
H <- 123
t$ring <- 100 * t$ring / H
t$troph <- 100 * (t$mid + t$troph) / H             # the 130 grey is the blend
t$schiz <- 100 * t$schiz / H
p <- subset(t, x >= 5 & x <= 318)                  # drop the y-axis columns
x_new <- min(p$x[p$x > 200 & p$ring > 5])          # next cycle's rings reappear
x0 <- 2; idc <- x_new - x0
frac <- function(xx) (xx - x0) / idc
lin <- function(xa, ya, xb, yb, y) xa + (y - ya) * (xb - xa) / (yb - ya)
## linear interpolation to the crossing, rather than the first column past it
i50 <- which(p$ring < 50)[1]
x50 <- lin(p$x[i50 - 1], p$ring[i50 - 1], p$x[i50], p$ring[i50], 50)
irt <- which(p$ring < p$troph)[1]
dprev <- p$ring[irt - 1] - p$troph[irt - 1]; dthis <- p$ring[irt] - p$troph[irt]
xrt <- lin(p$x[irt - 1], dprev, p$x[irt], dthis, 0)
i2 <- which(p$x > 50 & p$ring < 2)[1]
x2 <- lin(p$x[i2 - 1], p$ring[i2 - 1], p$x[i2], p$ring[i2], 2)

cat("\n\n2. LLINAS ET AL. (2006) Figure 1, 3D7 panel, read off the figure.\n",
    "   Peak ring ", sprintf("%.0f", max(p$ring)), "%; next cycle's rings reappear at ",
    sprintf("%.0f%%", 100 * frac(x_new)), " of the panel, taken as the IDC.\n\n", sep = "")
out <- data.frame(
    boundary = c("ring % falls below 50", "ring % crosses trophozoite %",
                 "ring % essentially 0 (<2)"),
    frac_of_IDC = round(c(frac(x50), frac(xrt), frac(x2)), 3),
    hours_at_48 = round(48 * c(frac(x50), frac(xrt), frac(x2)), 1))
print(out, row.names = FALSE)

cat("\n   Staging profile every 5% of the IDC (percent of parasites):\n\n")
g <- data.frame(pct_IDC = seq(0, 100, by = 5))
g$x <- x0 + g$pct_IDC / 100 * idc
g$h48 <- round(48 * g$pct_IDC / 100, 1)
for (v in c("ring", "troph", "schiz"))
    g[[v]] <- round(approx(p$x, p[[v]], xout = g$x, rule = 2)$y)
print(g[, c("pct_IDC", "h48", "ring", "troph", "schiz")], row.names = FALSE)

cat("\n\n3. THE WINDOW FOR THE STARTING AGE DISTRIBUTION.\n",
    "   The two sources agree to within ~1 h on a 48 h cycle. At the cycle\n",
    "   lengths these data fit, the ring window is:\n\n", sep = "")
for (cl in c(39.4, 41.0, 42.844, 45.6))
    cat(sprintf("   cycle_length %5.2f h: Smith 0-%.1f h, Llinas (ring < 2%%) 0-%.1f h\n",
                cl, 0.312 * cl, frac(x2) * cl))
cat("\n   Against the CIRCULATING window the design currently uses (24 h at a\n",
    "   48 h cycle; half-sequestration 18.58 h, 5% circulating 24.3 h), the\n",
    "   ring window is about 60-65% as wide.\n", sep = "")
