using DSGE

m = AnSchorfheide()

df = load_data(m; verbose=:high)

estimate(m, df; verbose=:high)
