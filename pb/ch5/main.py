### Eli Blaney ###
###  BIO 501   ###
### Chapter 5  ###
### 2022-01-24 ###
print("Chapter 5 Exercises\n")

# Regular expressions package for my satisfaction
import re

### Exercise 5.1a ###
print("\nExercise 5.1a")

def clean_sequence(sequence, validChars="GPAVLIMCFYWHKRQNEDST"):
    """Cleans sequences by removing non-whitelisted characters and converting to uppercase.

    Keyword arguments:
    sequence -- The sequence to clean.
    validChars -- Characters to keep. Default includes regular amino acid letter codes.

    """

    # To be case-insensitive, both the sequence and validChars must be transformed to uppercase.
    sequence = sequence.upper()
    validChars = validChars.upper()

    # Build the regular expressions pattern to exlude invalid characters
    # [^abc] will match all characters except "a", "b", or "c"
    # In a similar way, [^validChars] will match all characters except our whitelisted residues
    pattern = r"[^" + validChars + "]"

    # Now, we can use re.sub() to substitute all non-whitelisted characters with an empty string
    return re.sub(pattern, "", sequence)

def residue_percentage(sequence, residues="AILMFWYV", decimals=0):
    """Determines the percentage of residues in a sequence.

    Keyword arguments:
    sequence -- The larger sequence containing characters to count. Does not necessarily need to be a protein sequence.
    residues -- The residues to count in the larger sequence. Default includes hydrophobic residue letter codes.

    """

    # First clean the sequence by keeping only characters that are amino acid residue codes
    sequence = clean_sequence(sequence)
    # We can now safely use the length of this cleaned string in the calculations
    total = len(sequence)

    # If the string is completely invalid, we would iterate unnecessarily and also get a divide by 0 error
    # We can avoid that by having an early return statement if there are no valid characters in the sequence
    if total == 0:
        return 0

    # content is a counter variable for all the characters that we want to match, cumulatively
    content = 0
    # Iterate through all of the residues. Because of this iteration, residues can be a string or a list
    for r in residues:
        # Apply upper() to each character in the string (or to each string in the list)
        # Then, count its appearances in the cleaned sequence and add it to the content counter
        content += sequence.count(r.upper())

    # Divide the counted characters by the total length and convert to a percentage
    # Then round to a default of 0 decimal places and return
    return round(content / total * 100, decimals)

# Test basic functionality
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", "M") == 5
# Test lowercase character case
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", "r") == 10
# Test basic functionality again
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", "L") == 50
# Test missing residue case
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", "Y") == 0

print("Exercise 5.1a finished without errors.")

### Exercise 5.1b ###
print("\nExercise 5.1b")

# Test basic functionality
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", ["M"]) == 5
# Test multiple residues case
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", ['M', 'L']) == 55
# Test even more residues case
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP", ['F', 'S', 'L']) == 70
# Test default hydrophobic argument case
assert residue_percentage("MSRSLLLRFLLFLLLLPPLP") == 65

print("Exercise 5.1b finished without errors.")