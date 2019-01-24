DIR='../../data/fullPhysicsTest/run52/plots/plan/'
convert -loop 0 -delay 10 ${DIR}scatter_2d_XZ_*.png ${DIR}anim_xz.gif
convert -loop 0 -delay 10 ${DIR}scatter_2d_XY_*.png ${DIR}anim_xy.gif
convert -loop 0 -delay 10 ${DIR}scatter_3d_*.png ${DIR}anim_xyz.gif
