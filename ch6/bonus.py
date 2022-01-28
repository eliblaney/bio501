### Bonus Exercise 6.6 ###
print("\nBonus Exercise 6.6\n")

# Find the factorial using a for loop
def factorial_for(x):
    # For inputs outside the domain of this function, we can return None
    if x < 0:
        return None

    # This product variable will be repeatedly multiplied
    # by each required number to achieve the factorial.
    product = 1
    # Iterate through each number from 0 to x-1 (inclusive)
    for xi in range(x):
        # Set and multiply the product by each number from 1-x (inclusive)
        product *= (xi + 1)

    # We are left with the factorial result
    return product

# Find the factorial using a while loop
def factorial_while(x):
    # For inputs outside the domain of this function, we can return None
    if x < 0:
        return None

    # This counter variable keeps track of each number from 1 to x (inclusive)
    i = 1
    # Like in factorial_for(), this product variable will be multiplied to achieve
    # the factorial.
    product = 1
    # Iterate from 1 to x (inclusive) using a while loop
    while i <= x:
        # Multiply the product by each number in the sequence
        product *= i
        # Update the counter
        i += 1

    # We are left with the factorial
    return product

# Find the factorial using recursion
def factorial_recursion(x):
    # For inputs outside the domain of this function, we can return None
    if x < 0:
        return None

    # Base case to stop recursion upon reaching 0 or 1
    if x == 0 or x == 1:
        return 1

    # Multiply x by the number below x and recursively repeat this process
    # for every number down to 1 (or 0).
    return x * factorial_recursion(x - 1)

# Test it out with 5! = 120
x = 5
print("Factorial of", x, "using a for loop is", factorial_for(x))
print("Factorial of", x, "using a while loop is", factorial_while(x))
print("Factorial of", x, "using recursion is", factorial_recursion(x))
