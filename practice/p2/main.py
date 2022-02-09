# Iterate through all of the cell concentration test files.
for i in range(4):
    # Keep track of the filename for our printed output.
    filename = 'cell_Concentrations' + str(i + 1) + '.txt'
    file = open(filename)

    # These reasons are "English" and represent reasons for why
    # ENG101 would not have been created.
    reasons = []
    # Go through all of the given variables in the file.
    for line in file:
        # Split up the line into two variables that keep
        # track of the element (e.g. CSC121) and whether
        # it is present or not.
        element, is_present = line.rstrip().split(",")
        # For the inhibitors specifically, all we need to
        # check is if they are present, as that causes an issue.
        if element == "Na+" or element == "CSC121":
            if is_present == " True":
                # If they are, just say that they were present.
                reasons.append(element + " was present")
        # For everything else, we need to check if they are absent.
        elif is_present == " False":
            # If they are, make sure we detail that.
            reasons.append("no " + element + " present")

    # If there are any reasons that we couldn't make ENG101, then
    # we can print those out. Although reasons is a list, it is
    # "truthy" when it has elements inside it and "falsey" when
    # it is empty.
    if reasons:
        print("[" + filename + "] ENG101 could not be made because " + ", ".join(reasons))
    else:
        print("[" + filename + "] ENG101 was made.")

    # We don't need the file anymore; continue to the next one.
    file.close()