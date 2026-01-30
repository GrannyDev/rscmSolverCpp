import random

INPUT_FILE = "0&1&2_adders_constants.txt"
MASK_INPUT_FILE = "0&1_adders_constants.txt"
OUTPUT_FILE = "dataSets/2ADD/timeout/12bits/256on12bitsR100.txt"
NUM_LINES = 100
NUMBERS_PER_LINE = 256

def load_numbers(path):
    """Reads numbers from a file and filters those in [0, 31]."""
    numbers = set()
    with open(path, "r") as file:
        for line in file:
            num_str = line.split(":")[0].strip()
            if num_str.isdigit():
                num = int(num_str)
                if 0 <= num <= 4095:
                    numbers.add(num)
    return list(numbers)

def generate_random_lines(valid_numbers, mask_numbers):
    """
    Creates 100 lines of randomly selected numbers where each line contains at
    least one value not present in the mask file.
    """
    if len(valid_numbers) < NUMBERS_PER_LINE:
        raise ValueError("Not enough valid numbers in range 0-31.")

    non_mask_numbers = [n for n in valid_numbers if n not in mask_numbers]
    if not non_mask_numbers:
        raise ValueError("All valid numbers are masked; cannot satisfy constraint.")

    with open(OUTPUT_FILE, "w") as file:
        for _ in range(NUM_LINES):
            while True:
                selected_nums = random.sample(valid_numbers, NUMBERS_PER_LINE)
                if any(n not in mask_numbers for n in selected_nums):
                    break
            file.write(";".join(map(str, selected_nums)) + "\n")

def main():
    valid_numbers = load_numbers(INPUT_FILE)
    mask_numbers = set(load_numbers(MASK_INPUT_FILE))
    generate_random_lines(valid_numbers, mask_numbers)
    print(f"Generated {NUM_LINES} lines in {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
