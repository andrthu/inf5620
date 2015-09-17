from wave1D_solver import *
from numpy import *


def viz(
    I, V, f, q, L, dt, C, T, ue, max_q,  # PDE paramteres
    umin, umax,               # Interval for u in plots
    animate=True,             # Simulation with animation?
    tool='matplotlib',        # 'matplotlib' or 'scitools'
    solver_function=solver,   # Function with numerical algorithm
    plot_error=False,
    disc='type57'
    ):
    """Run solver and visualize u at each time level."""
    sleep=0.002
    def plot_u_st(u, x, t, n):
        """user_action function for solver."""
        plt.plot(x, u, 'r-',
                 xlabel='x', ylabel='u',
                 axis=[0, L, umin, umax],
                 title='t=%f' % t[n], show=True)
        # Let the initial condition stay on the screen for 2
        # seconds, else insert a pause of 0.2 s between each plot
        time.sleep(2) if t[n] == 0 else time.sleep(sleep)
        plt.savefig('frame_%04d.png' % n)  # for movie making

    class PlotMatplotlib:
        def __call__(self, u, x, t, n):
            """user_action function for solver."""
            if plot_error==True:
                if n == 0:
                    plt.ion()
                    self.lines = plt.plot(x, u-ue(x,0), 'r-')
                    plt.xlabel('x');  plt.ylabel('error')
                    plt.axis([0, L, umin, umax])
                    plt.legend(['t=%f' % t[n]], loc='lower left')
                else:
                    self.lines[0].set_ydata(u-ue(x,t[n]))
                    plt.legend(['t=%f' % t[n]], loc='lower left')
                    plt.draw()
                time.sleep(2) if t[n] == 0 else time.sleep(sleep)
                plt.savefig('tmp_%04d.png' % n)  # for movie making
            else:
                if n == 0:
                    plt.ion()
                    self.lines = plt.plot(x, u, 'r-')
                    self.lines2 = plt.plot(x,ue(x,0),'b--')
                    plt.xlabel('x');  plt.ylabel('error')
                    plt.axis([0, L, umin, umax])
                    plt.legend(['t=%f' % t[n]], loc='lower left')
                else:
                    self.lines[0].set_ydata(u)
                    self.lines2[0].set_ydata(ue(x,t[n]))
                    plt.legend(['t=%f' % t[n]], loc='lower left')
                    plt.draw()
                    
                time.sleep(2) if t[n] == 0 else time.sleep(sleep)
                plt.savefig('tmp_%04d.png' % n)  # for movie making




    if tool == 'matplotlib':
        import matplotlib.pyplot as plt
        plot_u = PlotMatplotlib()
    elif tool == 'scitools':
        import scitools.std as plt  # scitools.easyviz interface
        plot_u = plot_u_st
    import time, glob, os

    # Clean up old movie frames
    for filename in glob.glob('tmp_*.png'):
        os.remove(filename)

    # Call solver and do the simulaton
    user_action = plot_u if animate else None
    u, x, t, cpu = solver_function(
        I, V, f, q, L, dt, C, T, max_q,  user_action, disc)

    # Make video files
    fps = 4  # frames per second
    codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                     libtheora='ogg')  # video formats
    filespec = 'tmp_%04d.png'
    movie_program = 'ffmpeg'  # or 'avconv'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s movie.%(ext)s' % vars()
        os.system(cmd)

    if tool == 'scitools':
        # Make an HTML play for showing the animation in a browser
        plt.movie('tmp_*.png', encoder='html', fps=fps,
                  output_file='movie.html')
    return cpu


if __name__ == "__main__":

    ue=lambda x,t: sin(pi*x)*cos(pi*t)

    viz(lambda x:sin(pi*x),None,None, 1,1,0.1,0.9,20,
        ue,-2,2,True,'matplotlib',solver,True)


"""
    I, V, f, c, L, dt, C, T,ue,  
    umin, umax,               
    animate=True,             
    tool='matplotlib',        
    solver_function=solver,   
    )
"""
