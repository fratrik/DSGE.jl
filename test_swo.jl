using DSGE

m = SmetsWoutersOrig()

df = load_data(m; verbose=:high)

estimate(m, df; verbose=:high)
