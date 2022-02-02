# Open the file with the book entries
books = open('BookInfo.txt')

# This list will hold the books that we want to keep
filtered_books = []
# Iterate through each line in the file
for entry in books:
    # Split the line up into three distinct variables depending on what they represent
    title, author, pages = entry.split(", ")
    # Check if it's a book that we want
    if title != "The First World War" and int(pages) > 500 and int(pages) < 800:
        # Add the books we like to the list and remove the newline character
        filtered_books.append(entry.rstrip())

# Sort the list by author's last name, making use of a lambda function.
# The lambda works by telling Python that to find the item to sort the list
# by, it should split each entry up with ", " and look at the second element (index 1).
# That element is the author's last name, which will be automatically sorted alphabetically.
filtered_books.sort(key=lambda entry: entry.split(", ")[1])

# Print out our list! Problem requests the entries, so they are printed as-is.
print(filtered_books)
