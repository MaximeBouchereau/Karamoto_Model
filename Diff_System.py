import numpy as np
import matplotlib.pyplot as plt
import os


class Integration:
    """Class for numerical integration with Kuramoto model."""

    @staticmethod
    def F(Theta, K, omega):
        """Vector field associated to Kuramoto model.
        Input:
        - Theta: Aray of shape (d,1) - d is the dimensionality of the problem - Space variable.
        - K: Float - Coupling parameter.
        - omega: Array of shape (d,) - Proper frequencies of each angle."""

        Theta = Theta.reshape(Theta.shape[0], 1)
        d = Theta.shape[0]
        V = Theta @ np.ones((1, d)) - np.ones((d, 1)) @ Theta.T

        FF = omega - K * np.mean(np.sin(V), axis=1)

        return FF

    @staticmethod
    def RP(Theta):
        """Similarity between angles. Near 0: differences. Near 1: similarities.
        Input:
        - Theta: Aray of shape (d,) - d is the dimensionality of the problem - Space variable."""

        R = (np.mean(np.cos(Theta)) ** 2 + np.mean(np.sin(Theta)) ** 2) ** 0.5
        return R

    @staticmethod
    def Integrate(params, save = False):
        """Integrate Kuramoto model.
        Inputs:
        - params: Dict - Dictionary of parameter values:

            > name: Str - Name of model.
            > d: Int - Dimensionality of the problem.
            > K: Float - Coupling parameter.
            > mu: Float - Expected value for initial distribution of angular velocities distribution
            > sigma: Float - Standard deviation for initial distribution of angular velocities distribution
            > T: Float - Length of time interval resolution ([0, T]).
            > delta_t: Float - Step size.
            > num_meth: Str - Numerical method for integration. "FEuler" or "RK2".

        - save: Bool - Saves the figure or not. Default: False.

        => Angular velocities are distributed with normal distribution.
        => Distribution of initial conditions is uniform of [0, 2+pi].
        """

        name, d, K, mu, sigma, T, delta_t, num_meth = params['name'], params['d'], params['K'], params['mu'], params['sigma'], params['T'], params['delta_t'], params['num_meth']

        if num_meth not in ["FEuler", "RK2"]:
            raise ValueError("num_meth: must be either FEuler or RK2")

        omega = np.random.normal(1, 0.01, (d,))

        TT = np.arange(0, T + delta_t, delta_t)
        SYNCH = np.zeros_like(TT)
        N = len(TT)
        THETA = np.zeros((d, N))

        # INITIALIZATION
        THETA[:, 0] = np.random.uniform(0, 2 * np.pi, size=(d,))
        SYNCH[0] = Integration.RP(THETA[:, 0])

        # NUMERICAL INTEGRATION
        for n in range(N-1):
            print(" > Resolution: " + str(n+1) + " / " + str(N-1) + "  -  R = " + format(Integration.RP(THETA[:, n]), '.3f'), end="\r")
            if num_meth == "FEuler":
                THETA[:, n + 1] = THETA[:, n] + delta_t * Integration.F(THETA[:, n], K, omega)
            if num_meth == "RK2":
                THETA[:, n + 1] = THETA[:, n] + delta_t * Integration.F(THETA[:, n] + delta_t / 2 * Integration.F(THETA[:, n], K, omega), K, omega)
            SYNCH[n + 1] = Integration.RP(THETA[:, n + 1])


        # PLOT OR SAVE RESULTS
        if save == False:
            plt.figure()
            plt.plot(TT, np.mod(THETA.T, 2 * np.pi))
            plt.xlabel("$t$")
            plt.ylabel(r"$\theta$")
            plt.title(r"$\theta$")
            plt.grid()

            plt.figure()
            plt.plot(TT, SYNCH)
            plt.xlabel("$t$")
            plt.ylabel("$R(t)$")
            plt.title("Angle similarity")
            plt.grid()
            plt.show()

        if save == True:
            name_file = params['name'] + "_"
            name_file += "d=" + str(params['d']) + "_"
            name_file += "K=" + str(params['K']) + "_"
            name_file += "mu=" + str(params['mu']) + "_"
            name_file += "sigma=" + str(params['sigma']) + "_"
            name_file += "T=" + str(params['T']) + "_"
            name_file += "DeltaT=" + str(params['delta_t']) + "_"
            name_file += "NumMeth=" + str(params['num_meth'])
            os.makedirs(name_file, exist_ok=False)
            np.savez(name_file + "/Integration.npz", solution=THETA, direction_match=SYNCH, parameters=params)

        return None
