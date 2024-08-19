import numpy as np

# The original number
number = 9.5367431640625e-7

# Square the number
squared_number = number ** 2

# Convert the squared number to its bit representation using numpy
bit_representation = np.binary_repr(np.float64(squared_number).view(np.int64), width=64)
print(bit_representation)

# Extract the original bit representation components
sign_bit = bit_representation[0]
exponent_bits = bit_representation[1:12]
mantissa_bits = bit_representation[12:]

# Convert the exponent to an integer
exponent_value = int(exponent_bits, 2)

# Increase the exponent by 4
new_exponent_value = exponent_value - 4

# Ensure the new exponent is still within the range of 11 bits
new_exponent_bits = format(new_exponent_value, '011b')

# Construct the new bit representation
new_bit_representation = sign_bit + new_exponent_bits + mantissa_bits
print(new_bit_representation)


# Convert the new bit representation back to a floating-point number
new_number = np.int64(int(new_bit_representation, 2)).view(np.float64)

# Calculate the square root of the new number
square_root = np.sqrt(new_number)

print(new_number, square_root)
