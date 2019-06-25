#DIR='../../data/fullPhysicsTest/run51/plots/plan/anim1/'
#convert -loop 0 -delay 10 ${DIR}scatter_2d_XZ_*.png ${DIR}anim_xz.gif
#convert -loop 0 -delay 10 ${DIR}scatter_2d_XY_*.png ${DIR}anim_xy.gif
#convert -loop 0 -delay 10 ${DIR}scatter_3d_*.png ${DIR}anim_xyz.gif

#DIR='../../data/fullPhysicsTest/run51/plots/plan/anim1/'
#convert -loop 0 -delay 10 ${DIR}anim_*.png ${DIR}1anim.gif


DIR='../../data/prodRuns/run100/plots/planAnim_scatter/'
convert -loop 0 -delay 25 ${DIR}anim_*.png ${DIR}anim.gif


#DIR='../../data/prodRuns/run100/plots/planAnim_scatter'
#ffmpeg -framerate 5 -i $DIR/anim_%03d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $DIR/anim.mp4
#ffmpeg -framerate 5 -i $DIR/anim_*.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $DIR/anim.mp4
