#!/usr/bin/env python3
import os
import fitz  # PyMuPDF
import re
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

def clean_text(text):
    """Clean and format extracted text"""
    # Remove excessive whitespace
    text = re.sub(r'\s+', ' ', text)
    # Fix common PDF extraction issues
    text = re.sub(r'(?<=[a-z])(?=[A-Z])', ' ', text)  # Add space between camelCase
    # Preserve newlines for paragraphs
    text = re.sub(r'\s*\n\s*\n\s*', '\n\n', text)
    return text.strip()

def extract_images(page, page_num, output_dir):
    """Extract images from a PDF page"""
    image_list = page.get_images()
    images = []
    
    for img_index, img in enumerate(image_list):
        # Get the image
        xref = img[0]
        pix = fitz.Pixmap(page.parent, xref)
        
        # Convert to RGB if necessary
        if pix.n - pix.alpha < 4:
            pix = fitz.Pixmap(fitz.csRGB, pix)
        
        # Save image
        img_filename = f"page{page_num}_img{img_index}.png"
        img_path = os.path.join(output_dir, img_filename)
        pix.save(img_path)
        images.append(img_filename)
        pix = None
    
    return images

def pdf_to_markdown(pdf_path, output_dir):
    """Convert a single PDF to markdown"""
    pdf_name = Path(pdf_path).stem
    lecture_dir = Path(pdf_path).parent.name
    
    # Create output directory for this PDF
    pdf_output_dir = os.path.join(output_dir, lecture_dir, pdf_name)
    os.makedirs(pdf_output_dir, exist_ok=True)
    
    # Create images directory
    images_dir = os.path.join(pdf_output_dir, "images")
    os.makedirs(images_dir, exist_ok=True)
    
    markdown_content = f"# {pdf_name}\n\n"
    
    try:
        # Open PDF
        doc = fitz.open(pdf_path)
        
        print(f"Processing {pdf_name} ({doc.page_count} pages)...")
        
        for page_num in range(doc.page_count):
            page = doc[page_num]
            
            # Add page header
            markdown_content += f"\n## Page {page_num + 1}\n\n"
            
            # Extract text
            text = page.get_text()
            if text.strip():
                cleaned_text = clean_text(text)
                markdown_content += cleaned_text + "\n\n"
            
            # Extract images
            images = extract_images(page, page_num + 1, images_dir)
            for img in images:
                markdown_content += f"![{img}](images/{img})\n\n"
            
            # Extract tables if any
            tables = page.find_tables()
            if tables:
                markdown_content += "### Tables\n\n"
                for table_idx, table in enumerate(tables):
                    markdown_content += f"Table {table_idx + 1}:\n\n"
                    # Convert table to markdown
                    for row in table.extract():
                        markdown_content += "| " + " | ".join(str(cell) if cell else "" for cell in row) + " |\n"
                    markdown_content += "\n"
        
        # Save markdown file
        md_path = os.path.join(pdf_output_dir, f"{pdf_name}.md")
        with open(md_path, 'w', encoding='utf-8') as f:
            f.write(markdown_content)
        
        doc.close()
        
        print(f"✓ Converted {pdf_name} → {md_path}")
        return True, pdf_name
        
    except Exception as e:
        print(f"✗ Failed to convert {pdf_name}: {str(e)}")
        return False, pdf_name

def main():
    # Find all PDFs
    pdf_files = []
    for root, dirs, files in os.walk("pdfs"):
        for file in files:
            if file.endswith('.pdf'):
                pdf_files.append(os.path.join(root, file))
    
    print(f"Found {len(pdf_files)} PDFs to convert\n")
    
    # Create output directory
    output_dir = "markdown_output"
    os.makedirs(output_dir, exist_ok=True)
    
    start_time = time.time()
    success_count = 0
    failed_files = []
    
    # Process PDFs with ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=4) as executor:
        # Submit all conversion tasks
        future_to_pdf = {
            executor.submit(pdf_to_markdown, pdf_path, output_dir): pdf_path 
            for pdf_path in pdf_files
        }
        
        # Process completed conversions
        for future in as_completed(future_to_pdf):
            success, filename = future.result()
            if success:
                success_count += 1
            else:
                pdf_path = future_to_pdf[future]
                failed_files.append(pdf_path)
    
    elapsed_time = time.time() - start_time
    
    print(f"\n{'='*50}")
    print(f"Conversion completed in {elapsed_time:.1f} seconds")
    print(f"Successfully converted: {success_count}/{len(pdf_files)} files")
    
    if failed_files:
        print(f"\nFailed conversions ({len(failed_files)}):")
        for pdf_path in failed_files:
            print(f"  - {pdf_path}")
    
    print(f"{'='*50}")

if __name__ == "__main__":
    main()