import fitz # PyMuPDF

doc = fitz.open("paper2.pdf")
print(f"Number of pages: {len(doc)}")

for i in range(len(doc)):
    page = doc[i]
    text = page.get_text()
    if "Spacer-Naive" in text or "degradation" in text:
        print(f"--- Page {i+1} ---")
        lines = text.split("\n")
        for line in lines:
            if "degradation" in line or "Spacer-Naive" in line or "evasion" in line or "Cascade" in line:
                print(line)
