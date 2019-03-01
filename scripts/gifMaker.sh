#DIR='../../data/fullPhysicsTest/run51/plots/plan/anim1/'
#convert -loop 0 -delay 10 ${DIR}scatter_2d_XZ_*.png ${DIR}anim_xz.gif
#convert -loop 0 -delay 10 ${DIR}scatter_2d_XY_*.png ${DIR}anim_xy.gif
#convert -loop 0 -delay 10 ${DIR}scatter_3d_*.png ${DIR}anim_xyz.gif

#DIR='../../data/fullPhysicsTest/run51/plots/plan/anim1/'
#convert -loop 0 -delay 10 ${DIR}anim_*.png ${DIR}1anim.gif


DIR='../../data/fullPhysicsTest/run51_hc/plots/masterAnim'
ffmpeg -framerate 10 -i $DIR/anim_%000d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $DIR/anim.mp4
