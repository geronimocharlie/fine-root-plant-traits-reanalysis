# Synthetic dataset to illustrate varimax rotation (PCA)
# - Creates two latent gradients (F1, F2)
# - Traits are mixtures of both -> unrotated PCs have mixed loadings
# - Varimax rotation makes loadings more "one-hot"
# ======================================================

set.seed(42)

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---------------------------
# 1) Build synthetic data X
# ---------------------------
n <- 180  # "cities" (points)
# latent factors (gradients)
F1 <- rnorm(n)
F2 <- rnorm(n)

# Create 6 traits as linear mixtures of F1 and F2 + noise
# First 3 traits lean toward F1, last 3 toward F2, but all are mixed
noise <- function(sd=0.35) rnorm(n, sd=sd)

T1 <- 0.90*F1 + 0.50*F2 + noise()
T2 <- 0.80*F1 + 0.40*F2 + noise()
T3 <- 0.70*F1 + 0.60*F2 + noise()

T4 <- 0.50*F1 + 0.90*F2 + noise()
T5 <- 0.40*F1 + 0.80*F2 + noise()
T6 <- 0.60*F1 + 0.70*F2 + noise()

X <- data.frame(T1,T2,T3,T4,T5,T6)

# Standardize like typical PCA preprocessing
Xz <- scale(X)

# ---------------------------
# 2) PCA + Varimax rotation
# ---------------------------
pca <- prcomp(Xz, center=FALSE, scale.=FALSE)

k <- 2
L  <- pca$rotation[, 1:k, drop=FALSE]  # loadings (traits x 2)
S  <- pca$x[, 1:k, drop=FALSE]         # scores   (points x 2)

# Varimax rotation on loadings
vr    <- varimax(L)
L_rot <- vr$loadings[,]
Rmat  <- vr$rotmat

# Rotate scores consistently
S_rot <- S %*% Rmat

# ---------------------------
# 3) Plot A: point cloud + two axis systems
# ---------------------------
# ---- PLOT A: use rotated scores so the cloud aligns with rotated axes ----
# Plot points in UNROTATED PCA score space
df_scores <- data.frame(PC1 = S[,1], PC2 = S[,2])

# Rotated axes directions in the unrotated coordinate system:
# columns of Rmat give where rotated axes point to (in PC coordinates)
ang_rot1 <- atan2(Rmat[2,1], Rmat[1,1])
ang_rot2 <- atan2(Rmat[2,2], Rmat[1,2])

Llen <- max(abs(c(df_scores$PC1, df_scores$PC2))) * 1.10

# helper: segments through origin
axis_segment <- function(angle_rad, length, type){
  data.frame(
    x  = -cos(angle_rad)*length,
    y  = -sin(angle_rad)*length,
    xend =  cos(angle_rad)*length,
    yend =  sin(angle_rad)*length,
    type = type
  )
}

# Unrotated axes are just x and y axes (PC1, PC2)
ax_unrot1 <- axis_segment(0,    Llen, "Unrotated")   # PC1 axis
ax_unrot2 <- axis_segment(pi/2, Llen, "Unrotated")   # PC2 axis

# Rotated axes are diagonals
ax_rot1 <- axis_segment(ang_rot1, Llen, "Rotated")
ax_rot2 <- axis_segment(ang_rot2, Llen, "Rotated")

ax_all <- rbind(ax_unrot1, ax_unrot2, ax_rot1, ax_rot2)

dark_green <- "#0B3D2E"


df_scores <- data.frame(x = S_rot[,1], y = S_rot[,2])  # <-- key change

# In the rotated coordinate system, the original PCA axes are given by t(Rmat)
# (because Rmat maps unrot -> rot; so inverse maps rot -> unrot)
Rinv <- t(Rmat)

ang_unrot1 <- atan2(Rinv[2,1], Rinv[1,1])
ang_unrot2 <- atan2(Rinv[2,2], Rinv[1,2])

axis_lines <- function(angle_rad, length=1) {
  data.frame(
    x = c(-cos(angle_rad)*length, cos(angle_rad)*length),
    y = c(-sin(angle_rad)*length, sin(angle_rad)*length)
  )
}

Llen <- max(abs(c(df_scores$x, df_scores$y))) * 1.10

# dashed = original PCA axes (in rotated coordinates)
ax_unrot1 <- axis_lines(ang_unrot1, Llen)
ax_unrot2 <- axis_lines(ang_unrot2, Llen)

# solid = rotated axes (now just x/y axes in rotated score space)
ax_rot1 <- data.frame(x=c(-Llen, Llen), y=c(0,0))
ax_rot2 <- data.frame(x=c(0,0), y=c(-Llen, Llen))

# Natural dark-green palette
dark_green <- "#0B3D2E"   # points + rotated axes
mid_green  <- "#2F6F4E"   # optional
light_paper <- "#F4F1E7"  # background-ish

p_axes <- ggplot(df_scores, aes(x, y)) +
  geom_point(alpha=0.7, size=1.8, color=dark_green) +
  # unrotated PCA axes (dashed)
  geom_line(data=ax_unrot1, aes(x=x,y=y), linetype="dashed", linewidth=0.9, color="grey35") +
  geom_line(data=ax_unrot2, aes(x=x,y=y), linetype="dashed", linewidth=0.9, color="grey35") +
  # varimax axes (solid)
  geom_line(data=ax_rot1, aes(x=x,y=y), linewidth=1.2, color=dark_green) +
  geom_line(data=ax_rot2, aes(x=x,y=y), linewidth=1.2, color=dark_green) +
  coord_equal() +
  theme_minimal(base_size=12) +
    labs(
    title="Synthetic example",
    subtitle="Solid = original PCA axes; Dashed = varimax-rotated axes",
    x="PC1", y="PC2"
  )



# ---- PLOT B: loadings heatmap using absolute values ----

df_load_abs <- df_load %>%
  mutate(Loading = abs(Loading))

p_load <- ggplot(df_load_abs, aes(Component, Trait, fill=Loading)) +
  geom_tile() +
  facet_wrap(~Set, nrow=1) +
  scale_y_discrete(limits=rev) +
  theme_minimal(base_size=12) +
  theme(panel.grid=element_blank()) +
  scale_fill_gradient(
    low  = "#F4F1E7",   # paper
    high = "#0B3D2E"    # dark natural green
  ) +
  labs(
    title="Absolute loadings: clearer 'one-hot' structure after varimax",
    subtitle="Color intensity = loading strength (|value|)",
    x=NULL, y=NULL, fill="|Loading|"
  )

(p_axes | p_load) +
  plot_layout(widths=c(1.1, 1.6)) +
  plot_annotation(title="Varimax rotation: geometry unchanged, interpretation improved")


# Save if you want
ggsave("Results/Figures/varimax_demo.png", width=12, height=5.5, dpi=300)
