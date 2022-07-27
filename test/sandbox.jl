# using Revise
# using Distributions
# using Gadfly
# using DPMNeal3

# # Common settings
# M = 50;
# N = 500;
# iter = 10000;
# warmup = 5000;

# # Example 1
# dy = Normal();
# y = rand(dy, N);
# m = DPMNeal3.DPMNormal(; y, M);
# ỹ, fchain = DPMNeal3.sample!(m; iter, warmup);
# plot(
#     layer(x = ỹ, y = pdf.(dy, ỹ), Geom.line, Theme(default_color = "blue")),
#     layer(x = ỹ, y = mean(fchain), Geom.line, Theme(default_color = "red"))
# )

# # Example 2
# dy = MixtureModel(Normal, [(2, 1.0), (-1.0, 0.25)], [0.5, 0.5]);
# y = rand(dy, N);
# m = DPMNeal3.DPMNormal(; y, M);
# ỹ, fchain = DPMNeal3.sample!(m; iter, warmup);
# plot(
#     layer(x = ỹ, y = pdf.(dy, ỹ), Geom.line, Theme(default_color = "blue")),
#     layer(x = ỹ, y = mean(fchain), Geom.line, Theme(default_color = "red"))
# )
