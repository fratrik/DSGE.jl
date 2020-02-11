using DSGE

m = CurdiaReis()
# m = AnSchorfheide()

df = load_data(m; verbose=:none)
