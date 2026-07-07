import re

with open("paper2.tex", "r") as f:
    content = f.read()

# 1. Replace Figure~\ref{label} with \figref{label}  (green in-manuscript figures)
content = re.sub(r'Figure~\\ref\{([^}]+)\}', r'\\figref{\1}', content)

# 2. Replace Supplementary Figure S## with \suppfig{S##}  (red supplementary figures)
content = re.sub(r'Supplementary Figure\s+(S\d+)', r'\\suppfig{\1}', content)

# 3. Replace Supplementary Table S## with \supptab{S##}  (red supplementary tables)
content = re.sub(r'Supplementary Table\s+(S\d+)', r'\\supptab{\1}', content)

# 4. Replace standalone "Supplementary Material" (in abstract/conclusions) with red color
content = content.replace(
    "Supplementary Material.",
    "\\textcolor{red}{Supplementary Material}."
)

with open("paper2.tex", "w") as f:
    f.write(content)

print("Done. Replacements made:")
# Count how many of each type
print(f"  \\figref: {len(re.findall(r'\\\\figref\\{', content))}")
print(f"  \\suppfig: {len(re.findall(r'\\\\suppfig\\{', content))}")
print(f"  \\supptab: {len(re.findall(r'\\\\supptab\\{', content))}")
