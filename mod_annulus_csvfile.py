'''
This script processes a CSV file containing annulus radii data, 
and replaces specific NaN values with 0.0 when followed by a value > 0. The modified data is then saved to the same file.

Steps:
1. Reads the input CSV file into a pandas DataFrame.
2. Iterates through each column to find the first valid (non-NaN) value.
3. If the first valid value is not in the first row, checks the preceding row for NaN.
4. Replaces the NaN in the preceding row with 0.0.
5. Writes the modified DataFrame back to the existing CSV file, preserving the original header and excluding the index.
'''

# import pandas as pd

# csv_filename = "annulus_radii.csv"
# # Read CSV; header row (0) becomes column names automatically
# df = pd.read_csv(csv_filename)

# # Process each column except the frequency column
# # for col in df.columns.drop("frequency[GHz]"):
# for col in df.columns:
#     first_valid = df[col].first_valid_index()
#     # only if there's a valid value and it's not in the very first data row
#     if first_valid is not None and first_valid > 0:
#         prev_row = first_valid - 1
#         # replace that single NaN with 0.0
#         if pd.isna(df.at[prev_row, col]):
#             df.at[prev_row, col] = 0.0

# # Write back the CSV with the original header intact, no extra index
# df.to_csv(csv_filename, index=False, na_rep="nan")


#================================================================

import pandas as pd

df = pd.read_csv("annulus_radii.csv")

# if *_start has data and the next column is NaN.
# set to arbitrary large value
for col in df.filter(like="_start"):
    next_col = df.columns[df.columns.get_loc(col) + 1]
    mask = df[col].notna() & df[next_col].isna()
    df.loc[mask, next_col] = 0.013

# if *_end has data and the previous column is NaN.
# set to zero (begin center of FPI)
for col in df.filter(like="_end"):
    prev_col = df.columns[df.columns.get_loc(col) - 1]
    mask = df[col].notna() & df[prev_col].isna()
    df.loc[mask, prev_col] = 0
df.to_csv("annulus_radii.csv", index=False, na_rep="nan")
