using CSV, DataFrames

# read data for all operations
# as error take the last measurement, as it accumulate all errors from previous steps
# as runtime take the mean value over all data
table = CSV.read("out/efficiency/addition.csv", DataFrame)
addition_cipher_error = table.error_ciph_by_ciph_uncorrelated[1]
addition_plain_error = table.error_ciph_by_plain[1]
addition_scalar_error = table.error_ciph_by_scalar[1]
N_time = Int(table.N_time[1])
addition_cipher_time = sum(table.time_ciph_by_ciph[1:N_time])/N_time
addition_plain_time = sum(table.time_ciph_by_plain[1:N_time])/N_time
addition_scalar_time = sum(table.time_ciph_by_scalar[1:N_time])/N_time

table = CSV.read("out/efficiency/bootstrapping.csv", DataFrame)
bootstrap_error = table.error[1]
N_time = Int(table.N_time[1])
bootstrap_time = sum(table.time[1:N_time])/N_time

table = CSV.read("out/efficiency/double_bootstrapping.csv", DataFrame)
double_bootstrap_error = table.error[1]
N_time = Int(table.N_time[1])
double_bootstrap_time = sum(table.time[1:N_time])/N_time

table = CSV.read("out/efficiency/encryption-decryption.csv", DataFrame)
encryption_decryption_error = table.error[1]
N_time = Int(table.N_time[1])
encrypt_time = sum(table.time_encrypt[1:N_time])/N_time
decrypt_time = sum(table.time_decrypt[1:N_time])/N_time

table = CSV.read("out/efficiency/multiplication.csv", DataFrame)
multiplication_cipher_error = table.error_ciph_by_ciph_uncorrelated[1]
multiplication_plain_error = table.error_ciph_by_plain[1]
multiplication_scalar_error = table.error_ciph_by_scalar[1]
N_time = Int(table.N_time[1])
multiplication_cipher_time = sum(table.time_ciph_by_ciph[1:N_time])/N_time
multiplication_plain_time = sum(table.time_ciph_by_plain[1:N_time])/N_time
multiplication_scalar_time = sum(table.time_ciph_by_scalar[1:N_time])/N_time

table = CSV.read("out/efficiency/rotate.csv", DataFrame)
rotate_error = (table.error_rotate_1[1] + table.error_rotate_5[1] + 
                table.error_rotate_25[1])/3
N_time = Int(table.N_time[1])
rotate_time = (sum(table.time_rotate_1[1:N_time]) + sum(table.time_rotate_5[1:N_time]) + 
               sum(table.time_rotate_25[1:N_time]))/(3*N_time)

# use ciphertext addition as unit
error_unit = addition_cipher_error
time_unit = addition_cipher_time

# write out
table = DataFrame(encrypt_error = encryption_decryption_error/error_unit,
                  encrypt_time = encrypt_time/time_unit,
                  decrypt_error = encryption_decryption_error/error_unit,
                  decrypt_time = decrypt_time/time_unit,
                  addition_cipher_error = addition_cipher_error/error_unit,
                  addition_cipher_time = addition_cipher_time/time_unit,
                  addition_plain_error = addition_plain_error/error_unit,
                  addition_plain_time = addition_plain_time/time_unit,
                  addition_scalar_error = addition_scalar_error/error_unit,
                  addition_scalar_time = addition_scalar_time/time_unit,
                  multiplication_cipher_error = multiplication_cipher_error/error_unit,
                  multiplication_cipher_time = multiplication_cipher_time/time_unit,
                  multiplication_plain_error = multiplication_plain_error/error_unit,
                  multiplication_plain_time = multiplication_plain_time/time_unit,
                  multiplication_scalar_error = multiplication_scalar_error/error_unit,
                  multiplication_scalar_time = multiplication_scalar_time/time_unit,
                  rotate_error = rotate_error/error_unit,
                  rotate_time = rotate_time/time_unit,
                  bootstrap_error = bootstrap_error/error_unit,
                  bootstrap_time = bootstrap_time/time_unit,
                  double_bootstrap_error = double_bootstrap_error/error_unit,
                  double_bootstrap_time = double_bootstrap_time/time_unit)

CSV.write("out/efficiency/final_table_data.csv", table)
