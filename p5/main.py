plants = open('plants.txt')

# Dictionaries for metrics we want to keep
# track of during our program.
water_needs_per_m2 = {}
plant_albedos = {}

# Iterate through all of the lines in our data file
for line in plants:
    # Separate data to separate variables
    plant, albedo, water_needs, density = line.rstrip().split(",")
    # water_needs is L, density is m^2 (not m^-2)
    # So we have to divide them
    water_needs_per_m2[plant] = float(water_needs) / float(density)
    # Store the albedo as well so that we have it for the last part
    plant_albedos[plant] = float(albedo)
    # Print the albedo data as requested
    print("{} has a reflectivity of {:.2f}".format(plant, float(albedo)))

# We don't need the file past this point.
plants.close()

# Print out our water needs dictionary as requested
print("\nWater needs per m2 for all plants:")
print(water_needs_per_m2)

# I chose these two plants for this last part
plant1 = 'Hardy begonia'
plant2 = 'Geranium'

# Calculate the average albedo and total water
# usage using our dictionaries
average_albedo = (plant_albedos[plant1] + plant_albedos[plant2]) / 2
# Total water usage for a 10m2 plot assumes each plant
# is taking up 5m2. So we can add and multiply because of factoring.
total_water_usage = 5 * (water_needs_per_m2[plant1]  + water_needs_per_m2[plant2])

# Print our final data. The .format() method is really handy if
# you decide to learn it. You can see that I use it all the time.
print("\n{} and {} have an average albedo of {:.2f} and a total water usage of {:.2f} L for a 10m2 plot.".format(plant1, plant2, average_albedo, total_water_usage))
