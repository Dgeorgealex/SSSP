from plot_time_results_csv import fit_power_law
import numpy as np

def main():
    x_vals = np.array([
        5 * 10**5,
        1 * 10**6,
        2 * 10**6,
        5 * 10**6,
        1 * 10**7,
        2 * 10**7,
    ])

    log4 = x_vals * np.log2(x_vals) ** 4

    fit = fit_power_law(x_vals, log4)
    print(fit)

if __name__ == "__main__":
    main()