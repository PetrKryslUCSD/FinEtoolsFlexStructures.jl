# T3FF
results = Any[1.18423e-03, 1.22990e-03, 1.26116e-03, 1.27733e-03, 1.28531e-03, 1.28934e-03, 1.29124e-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)
results = Any[4.87096e-03, 5.00927e-03, 5.13014e-03, 5.19262e-03, 5.22241e-03, 5.23738e-03, 5.24359e-03]
q1, q2, q3 = results[5:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)


# S4 4, 8, 16, 32, 64
results = [1.2950943E-03, 1.2951259E-03, 1.2951361E-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)


# S4R5 8, 16, 32, 64
results = [1.2952425E-03, 1.2951767E-03, 1.2951408E-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)

# S3 8, 16, 32,
results = [5.2430298E-03, 5.2481278E-03,  5.2494994E-03 ]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)

# STRI3 8, 16, 32, 64
results = [5.2451388E-03, 5.2486794E-03, 5.2497592E-03 ]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)


# STRI3 8, 16, 32, 
results = [1.2945211E-03, 1.2949878E-03, 1.2951204E-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)

# S3 8, 16, 32,
results = [1.2944770E-03, 1.2949688E-03, 1.2951002E-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)


# S4 8, 16, 32,
results = [5.2495436E-03, 5.2497736E-03, 5.2498655E-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)


# S4R5 8, 16, 32, 
results = [5.2511394E-03, 5.2502993E-03, 5.2498932E-03]
q1, q2, q3 = results[end-2:end]
@show (q2^2 - q1*q3) / (2*q2 - q1 - q3)
