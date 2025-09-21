import matplotlib.pyplot as plt

import src.solver.define_equilibria as define_equilibria



if __name__ == '__main__':
    # a, b, c = 1, 1, 0
    # psi, R = gs_equation.get_psi(1, 0.5, 100, a, c)
    # p = gs_equation.get_p(a, psi)
    # F = gs_equation.get_F(b, psi)
    # plt.plot(psi, F)
    # plt.savefig("test.png")
    # plt.close()

    define_equilibria.define_multiple_equilibria()