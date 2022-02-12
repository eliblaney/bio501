# It would be better to use a library like pandas or even a dictionary
# to represent grades across subjects for multiple students, but for
# the goals of this assignment I will use a 2D array.
grades = []

# The first task is to get all the subject names
subjects = []
# As a placeholder, we'll use the name of our Average column
subject = "Average"
# We know that the total length must be enough at least for "Average",
# so we'll set that as our preliminary max_len_subject variable.
max_len_subject = len(subject)
# Continue getting user input until they press enter with no input
while subject:
    subject = input("Enter subject name (or enter to finish): ")
    if subject:
        # If they submit something, add it to the list
        subjects.append(subject)
        # Keep track of the maximum length of the subjects
        length = len(subject)
        if length > max_len_subject:
            max_len_subject = length

# Add our average column at the very end
subjects.append("Average")

print("\n-----------------------------\n")

# Same thing as before, but now with students
students = []
student = "Average"
# We know that we'll at least want enough for the "Average" cell
max_len_student = len(student)
# Continue getting student names until they stop giving input
while student:
    student = input("Enter student name (or enter to finish): ")
    if student:
        # Add the student name and keep track of the max length
        students.append(student)
        length = len(student)
        if length > max_len_student:
            max_len_student = length

# Add the row for the bottom of the table
students.append("Average")

print("\n-----------------------------\n")

# Now we want to add all the student grades for each subject
# into our 2D array.
# We'll take advantage of the same for loop so that we can
# find the average grades the students, too.
# We'll splice using [:-1] to avoid the bottom "average" row.
for student in students[:-1]:
    student_grades = []
    # Keeping track of the average grade for each student
    average = 0
    # Again, skip the "Average" column, but get all the
    # grades for each student in each subject.
    for subject in subjects[:-1]:
        grade = float(input("Input grade for {} in {}: ".format(student, subject)))
        # Add the grade to our running average. We'll divide later
        average += grade
        # We want to transform the grade back into a string,
        # because it will make printing easier. We'll also
        # round to one decimal place as we do this process.
        student_grades.append("{:.1f}".format(grade))
    # Find the actual average (again, ignoring the "average" column)
    average /= len(subjects) - 1
    # Add the average grade string, rounded to one decimal
    student_grades.append("{:.1f}".format(average))
    # Add the row of grades and student average
    grades.append(student_grades)

# Now we want to calculate the subject averages.
# We start with an empty array, which we'll add to.
grades.append([])
# We go through all of the subject columns, including
# the "average" column, which we can treat as a "subject".
for i in range(len(subjects)):
    subject = subjects[i]
    average = 0
    # We want to find the grade that each student
    # obtained in this subject, so we iterate through
    # all the students, ignoring the last "average" row.
    for j in range(len(students) - 1):
        # Now we add the student's grade in this subject
        # to our average.
        average += float(grades[j][i])
    average /= len(students)
    # Append the final average string to the bottom
    # row of our table. -1 tells Python to look at
    # the final index.
    grades[-1].append("{:.1f}".format(average))

# Now we want to print our grades table.
# This can get confusing because displaying tabular
# data in the console without external libraries can
# look a little weird at first, so stick with me here.

def divider(vertical_edge=False):
    """Creates a divider

    Keyword arguments:
    vertical_edge -- If the divider is for a vertical edge (top or
                     bottom of the table). Defaults to False.
    """
    edge = "|"
    if vertical_edge:
        edge = "-"
    # We want enough dashes for our student name as well as
    # enough dashes to cover all of the subject columns.
    print(edge + "-" * ((max_len_subject + 3) * len(subjects) + max_len_student + 2) + edge)

# The rest of the code prints the 2D grades array.
divider(True)
# This is the top left cell, for the student name column.
# We could print "Name", but I like having a blank cell.
print("|".ljust(max_len_student + 3), end="")

# Now we print the header row.
for subject in subjects:
    # Print each subject, spacing things out to our max_len_subject
    # string.ljust() is a method that left-justifies a string by
    # padding it with spaces to the right side. In other words,
    # all the resulting print statements will print the exact
    # same number of characters, making tabular data output.
    print("| {} ".format(subject.ljust(max_len_subject)), end="")
# The right edge of our header row.
print("|")

# Now we print the rest of the data.
for i in range(len(grades)):
    divider()
    student_grades = grades[i]
    # Print out the student's name, making sure they're all the
    # same number of characters using ljust() with max_len_student
    print("| {} ".format(students[i].ljust(max_len_student)), end="")
    # Now we print all their grades, padding them so that they are
    # the same length as the column name.
    for j in range(len(student_grades)):
        print("| {} ".format(student_grades[j].ljust(max_len_subject)), end="")
    # The right edge of the row.
    print("|")
# Bottom divider
divider(True)