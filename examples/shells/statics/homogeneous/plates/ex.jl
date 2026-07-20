using FinEtools
[ Info: n = 16
# [ Info: w(center): -476.47358596
# [ Info: q1 Range: -2.9321191370969727e6 to 41946.044216159105
# [ Info: q2 Range: -6.64926231952028e6 to -9355.418950102288
# [ Info: n = 32
# [ Info: w(center): -476.96988523
# [ Info: q1 Range: -3.0087674119316586e6 to 56014.801923275205
# [ Info: q2 Range: -7.853538696251049e6 to -5179.683232591838
# [ Info: n = 64
# [ Info: w(center): -477.08494442
# [ Info: q1 Range: -3.069801108404807e6 to 49312.135963368855
# [ Info: q2 Range: -9.129796114841904e6 to -2552.9371476161627
# [ Info: n = 128
# [ Info: w(center): -477.11269297
# [ Info: q1 Range: -3.1519778940834184e6 to 56689.820982265876
# [ Info: q2 Range: -1.0404172616023982e7 to -1266.0645530657014
# [ Info: n = 256
# [ Info: w(center): -477.11949346
# [ Info: q1 Range: -3.1863638011683472e6 to 59434.43280143945
# [ Info: q2 Range: -1.1672441437659271e7 to -630.3779682147219

s1 = -3.069801108404807e6
s2 = -3.1519778940834184e6
s3 = -3.1863638011683472e6 

@show (s2^2-s1*s3)/(2*s2-s1-s3) # -3.2111047274883403e6


s1 = -9.129796114841904e6
s2 = -1.0404172616023982e7
s3 = -1.1672441437659271e7

@show (s2^2-s1*s3)/(2*s2-s1-s3) # no convergence

# test_convergence_quarter_thickness 
# [ Info: Support soft --------------------------------------------------
# [ Info: Simply supported square plate with uniform load, Q4RS, stab_alpha=0.1  
# [ Info: Mesh distortion: striped
# [ Info: thickness/length = 0.02
# [ Info: n = 16
# [ Info: w(center): -476.50793207
# [ Info: q1 Range: -2.8675473975441414e6 to 23120.37140612474
# [ Info: q2 Range: -5.389897992636184e6 to -10407.342492025085
# [ Info: n = 32
# [ Info: w(center): -476.9780065
# [ Info: q1 Range: -3.0085282553927107e6 to 35317.52601524203
# [ Info: q2 Range: -6.664536941883095e6 to -5076.051301980432
# [ Info: n = 64
# [ Info: w(center): -477.08692322
# [ Info: q1 Range: -3.1016977260360303e6 to 49156.766252022506
# [ Info: q2 Range: -7.944618765146187e6 to -2500.9544211751613
# [ Info: n = 128
# [ Info: w(center): -477.11315229
# [ Info: q1 Range: -3.155337006932666e6 to 55387.226501020996
# [ Info: q2 Range: -9.218170230319126e6 to -1239.9666778262863
# [ Info: n = 256
# [ Info: w(center): -477.11958947
# [ Info: q1 Range: -3.183928629478992e6 to 59427.20320534183
# [ Info: q2 Range: -1.0486238184453264e7 to -617.4952834398622


s1 = -7.944618765146187e6
s2 = -9.218170230319126e6
s3 = -1.0486238184453264e7

@show (s2^2-s1*s3)/(2*s2-s1-s3) # no convergence
