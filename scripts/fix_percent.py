import re

with open("paper2.tex", "r") as f:
    lines = f.readlines()

fixed_count = 0
result = []

# Lines that legitimately use % as LaTeX (commands or comments)
latex_ok_pattern = re.compile(
    r'\\AtBeginDocument\{%|\\end\{tabular\}%|\\resizebox\{.*\}\{!\}\{%'
)

for i, line in enumerate(lines):
    lineno = i + 1
    
    # Skip pure comment lines (lines starting with optional spaces then %)
    stripped = line.lstrip()
    if stripped.startswith('%'):
        result.append(line)
        continue
    
    # Skip lines with known LaTeX % usage
    if latex_ok_pattern.search(line):
        result.append(line)
        continue
    
    # Find all bare % not preceded by backslash, not inside math ($...$)
    # Strategy: replace unescaped % with \%  but not at end of line (comment position)
    # We process char by char
    new_line = []
    in_math = False
    j = 0
    changed = False
    
    while j < len(line):
        c = line[j]
        
        if c == '$':
            in_math = not in_math
            new_line.append(c)
        elif c == '\\':
            # Escaped char – consume next char too
            new_line.append(c)
            if j + 1 < len(line):
                j += 1
                new_line.append(line[j])
        elif c == '%':
            # This is a comment if nothing follows that's not just newline → check if it's mid-text
            # Any % that appears after real content is either intentional or a bug
            # Check if everything before this on the line (excluding leading space) is LaTeX content
            before = ''.join(new_line).rstrip()
            if before:  # There is text before → this % truncates the line → fix it
                new_line.append('\\%')
                changed = True
                fixed_count += 1
                # Skip rest of line content after bare % (it was a comment, now escaped)
                # But we should NOT skip rest — just fix the %
            else:
                new_line.append(c)  # It's a leading/pure comment % — keep as is
        else:
            new_line.append(c)
        j += 1
    
    if changed:
        print(f"Line {lineno}: FIXED")
        print(f"  Before: {line.rstrip()[:100]}")
        print(f"  After:  {''.join(new_line).rstrip()[:100]}")
    
    result.append(''.join(new_line))

with open("paper2.tex", "w") as f:
    f.writelines(result)

print(f"\nTotal fixes: {fixed_count}")
