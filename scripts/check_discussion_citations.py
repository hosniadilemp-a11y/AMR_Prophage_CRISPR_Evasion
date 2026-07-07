import fitz  # PyMuPDF
import re

doc = fitz.open("paper2.pdf")

# Pages 8-12 are the Discussion section based on compile output
print("Scanning Discussion section for in-text citations...\n")
citation_pattern = re.compile(r'\((\d+(?:,\s*\d+)*)\)')
all_citations = set()

for i in range(7, len(doc)):  # pages 8-12 (0-indexed = 7-11)
    page = doc[i]
    text = page.get_text()
    citations = citation_pattern.findall(text)
    if citations:
        print(f"=== Page {i+1} ===")
        for c in citations:
            nums = [int(n.strip()) for n in c.split(',')]
            all_citations.update(nums)
        print(f"  Found: {citations}")

print(f"\nAll citation numbers found in Discussion: {sorted(all_citations)}")
