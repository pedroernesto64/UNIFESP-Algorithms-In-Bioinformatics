with open("output1.txt", "w") as f:
    print("This line will go into the file.", file=f)
    print(f"You can also print variables, like a number: {42}", file=f)
    print("This line also goes to the file.", file=f)

