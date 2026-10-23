import numpy as np
import matplotlib.pyplot as plt

import sys
from pathlib import Path

class Print:
    """Class for print results"""

    @staticmethod
    def Plot_Result(sol_path, save = False):
        """Print results: angle positions and angle similarity.
        Inputs:
        - sol_path: Str - Path to printed solution.
        - save: Bool - Save the result. Default: False.
        """

        SOL = np.load(sol_path, allow_pickle=True)
        THETA = SOL['solution']
        SYNCH = SOL['direction_match']
        params = SOL['parameters'].item()

        TT = np.arange(0, params['T'] + params['delta_t'], params['delta_t'])

        if save == True:
            name_file = params['name'] + "_"
            name_file += "d=" + str(params['d']) + "_"
            name_file += "K=" + str(params['K']) + "_"
            name_file += "mu=" + str(params['mu']) + "_"
            name_file += "sigma=" + str(params['sigma']) + "_"
            name_file += "T=" + str(params['T']) + "_"
            name_file += "DeltaT=" + str(params['delta_t']) + "_"
            name_file += "NumMeth=" + str(params['num_meth'])

        plt.figure()
        plt.plot(TT, np.mod(THETA.T, 2 * np.pi))
        plt.xlabel("$t$")
        plt.ylabel(r"$\theta$")
        plt.title(r"$\theta$")
        plt.grid()

        if save == True:
            plt.savefig(name_file + "/name_file_THETA.pdf", dpi=1500)

        plt.figure()
        plt.plot(TT, SYNCH, color="forestgreen")
        plt.xlabel("$t$")
        plt.ylabel("$R(t)$")
        plt.title("Angle similarity")
        plt.grid()

        if save == True:
            plt.savefig(name_file + "/name_file_SYNCH.pdf", dpi=1500)

        if save == False:
            plt.show()


        plt.close()
        SOL.close()
        return None


    @staticmethod
    def Plot_Results_Evolution(sol_path, save = False):
        """Plots evolution of system [Metronome synchronization].
        Input:
        - sol_path: Str - Path to printed solution.
        - save: Bool - Save the result. Default: False.
        """

        SOL = np.load(sol_path, allow_pickle=True)
        THETA = SOL['solution']
        SYNCH = SOL['direction_match']
        params = SOL['parameters'].item()

        name_file = params['name'] + "_"
        name_file += "d=" + str(params['d']) + "_"
        name_file += "K=" + str(params['K']) + "_"
        name_file += "mu=" + str(params['mu']) + "_"
        name_file += "sigma=" + str(params['sigma']) + "_"
        name_file += "T=" + str(params['T']) + "_"
        name_file += "DeltaT=" + str(params['delta_t']) + "_"
        name_file += "NumMeth=" + str(params['num_meth'])

        TT = np.arange(0, params['T'] + params['delta_t'], params['delta_t'])
        N = len(TT)

        ### PLOT OF METRONOME SYNCHRONIZATION
        d = THETA.shape[0]
        lbda = 1
        while d > lbda * int(1.5 * lbda):
            lbda += 1
        lbda_bis = int(np.sqrt(2 * d / 3))
        m = int((lbda + 1) * int(1.5 * (lbda + 1)) / (lbda + 1))
        a = int(d / m)
        b = d % m
        xx = np.mod(np.linspace(0, d-1, d), m)
        yy = np.linspace(0, d-1, d) // m

        ### FIGURE
        dx = np.max(xx) - np.min(xx) + 2
        dy = np.max(-yy) - np.min(-yy) + 2
        data_ratio = dy / dx

        # plt.figure(figsize=(5, 5 * data_ratio))
        # fig, ax = plt.subplots(figsize=(10, 10 * data_ratio))
        fig, ax = plt.subplots()

        plt.axis("equal")

        for n in range(N):
            print(" > Print: " + str(n) + " / " + str(N - 1), end="\r")
            if n % 1 == 0:
                ax.cla()
                ax.axis("off")


                ax.set_aspect('equal', adjustable='box')
                ax.set_xlim(np.min(xx) - 1, np.max(xx) + 1)
                ax.set_ylim(np.min(-yy) - 1, np.max(-yy) + 1)



                xx_b = xx + 0.5 * np.sin(np.cos(10 * THETA[:, n]))
                yy_b = yy - 0.5 * np.cos(np.cos(10 * THETA[:, n]))

                x_plot = np.vstack([xx, xx_b, np.full(d, np.nan)]).T.flatten()
                y_plot = np.vstack([yy, yy_b, np.full(d, np.nan)]).T.flatten()


                ax.plot(x_plot, -y_plot, color='blue', linestyle='-', label='Segments')
                ax.scatter(xx, -yy, color='red', zorder=3, label='(xx, yy)')
                ax.scatter(xx_b, -yy_b, color='green', zorder=3, label='(xx_b, yy_b)')

                # ax.scatter(xx + xx_b, - yy, color="r")
                # plt.xlim(-1, np.max(xx) + 1)
                # plt.ylim(np.min(-yy) - 1, np.max(-yy) + 1)
                ax.set_title(r"$t = $" + str(format(TT[n], '.2f')) + " - R(t) = " + str(format(SYNCH[n], '.2f')))

                evol_folder = Path(sol_path).parent / "Evolution"
                evol_folder.mkdir(parents=True, exist_ok=True)

                if save == True:
                    plt.draw()
                    plt.savefig(str(Path(sol_path).parent) + "/Evolution/evolution_" + str(n) + ".pdf", dpi=500)
                    plt.savefig(str(Path(sol_path).parent) + "/Evolution/evolution_" + str(n) + ".png", dpi=250)
                if save == False:
                    plt.pause(0.01)
                # plt.show()



        return None

