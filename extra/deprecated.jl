# function update_f!(m::AbstractDPM)
#     (; N, M, α, f, P, A, n) = skeleton(m)
#     α0, k0 = α[], first(P)
#     for i in 1:M
#         f[i] = α0 * exp(out_of_sample_logpredlik(m, i, k0))
#         for k in A 
#             f[i] += n[k] * exp(out_of_sample_logpredlik(m, i, k))
#         end
#         f[i] /= (N + α0)
#     end
#     return nothing
# end
