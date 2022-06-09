
import numpy as np
from ezc3d import c3d
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import seaborn as sns


trial_name = 'SautHaut_03'
# trial_name = '4TiretHaut_01'

c3d_experimental = c3d(trial_name+'.c3d')

#boucle pour supprimer les points non labelled
labels=c3d_experimental['parameters']['POINT']['LABELS']['value']
indices_supp=[]
for i in range (len(labels)) :
    if '*' in labels[i] :
        indices_supp=np.append(indices_supp,i)
ind_stop=int(indices_supp[0])

labels=c3d_experimental['parameters']['POINT']['LABELS']['value'][0:ind_stop]
# c3d_experimental['parameters']['POINT']['LABELS']['value']=labels




def plot_coutours():
    # Achanger pour ce qu'on veut vraiment comme frame
    ax.plot(np.array([0, 0.1]),
            np.array([0, 0.1]),
            np.array([0, 0.1]), '-w')
    return


def update(i, c3d_experimental, markers_point):
    for i_point in range(len(markers_point)):
        markers_point[i_point][0].set_data(np.array([c3d_experimental[0, i_point, i]]), np.array([c3d_experimental[1, i_point, i]]))
        markers_point[i_point][0].set_3d_properties(np.array([c3d_experimental[2, i_point, i]]))
    return


frame_range = [0, 9000]
# frame_range = [7000, 7020]
output_file_name = trial_name+'points.mp4'


fig_1 = plt.figure()
ax = p3.Axes3D(fig_1,auto_add_to_figure=False)
fig_1.add_axes(ax)
# ax = p3.Axes3D(fig_1)
ax.set_box_aspect([1, 1, 1])
# plot_coutours()

ax.set_xlim3d(
    [np.nanmin(c3d_experimental['data']['points'][0, :, :]), np.nanmax(c3d_experimental['data']['points'][0, :, :])])
ax.set_ylim3d(
    [np.nanmin(c3d_experimental['data']['points'][1, :, :]), np.nanmax(c3d_experimental['data']['points'][1, :, :])])
ax.set_zlim3d([0, np.nanmax(c3d_experimental['data']['points'][2, :, :])])

# colors_colormap = sns.color_palette(palette="viridis", n_colors=c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]].shape[1])
colors_colormap = sns.color_palette(palette="viridis", n_colors=c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]].shape[1])

# colors = [[] for i in range(c3d_experimental['data']['points'][0:ind_stop, :, frame_range[0]:frame_range[1]].shape[1])]
colors = [[] for i in range(c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]].shape[1])]
for i in range(c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]].shape[1]):
# for i in range(ind_stop):
    col_0 = colors_colormap[i][0]
    col_1 = colors_colormap[i][1]
    col_2 = colors_colormap[i][2]
    colors[i] = (col_0, col_1, col_2)

# markers_point = [ax.plot(0, 0, 0, '.', color=colors[i]) for i in range(c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]].shape[1])]
markers_point = [ax.plot(0, 0, 0, '.', color=colors[i]) for i in range(c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]].shape[1])]

# anim = animation.FuncAnimation(fig_1, update, frames=frame_range[1]-frame_range[0], fargs=(c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]], markers_point), blit=False)
anim = animation.FuncAnimation(fig_1, update, frames=frame_range[1]-frame_range[0], fargs=(c3d_experimental['data']['points'][:, 0:ind_stop, frame_range[0]:frame_range[1]], markers_point), blit=False)


anim.save(output_file_name, fps=20, extra_args=['-vcodec', 'libx264'])

plt.show()

