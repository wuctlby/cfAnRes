# pip install pypdf pymupdf
from pypdf import PdfReader, PdfWriter, Transformation
import fitz  # PyMuPDF
import math
from pathlib import Path

def collage_pdf_pages_to_single(
    input_pdf: str | list[str],
    output_pdf: str,
    output_png: str | None = None,
    page_width: float = 1920,      # 16:9 画布宽（point，1 point=1/72 inch）
    page_height: float = 1080,     # 16:9 画布高
    cols: int | None = None,       # 固定列数；不设则自动按 16:9 估计
    margin: float = 2.0,           # 每个小格四周留白；设更小=图间距更小
    render_zoom: float = 2.0       # 转 PNG 时的缩放（2.0 ≈ 144 DPI 基础上再放大）
):
    if isinstance(input_pdf, list):
        # 多个 PDF 合并为一个临时 PDF
        temp_writer = PdfWriter()
        for pdf_path in input_pdf:
            temp_reader = PdfReader(pdf_path)
            for page in temp_reader.pages:
                temp_writer.add_page(page)
        temp_path = Path(output_pdf).with_suffix('.temp.pdf')
        with open(temp_path, "wb") as f:
            temp_writer.write(f)
        input_pdf = str(temp_path)
    reader = PdfReader(input_pdf)
    n = len(reader.pages)
    if n == 0:
        raise ValueError("输入 PDF 没有页面。")

    # 自动估列数，尽量匹配 16:9
    aspect = page_width / page_height
    if cols is None:
        cols = math.ceil(math.sqrt(n * aspect))
    rows = math.ceil(n / cols)

    cell_w = page_width / cols
    cell_h = page_height / rows

    writer = PdfWriter()
    big = writer.add_blank_page(width=page_width, height=page_height)

    for i, src in enumerate(reader.pages):
        r, c = divmod(i, cols)

        # 源页尺寸
        sw = float(src.mediabox.width)
        sh = float(src.mediabox.height)

        # 可用区域（留出 margin 作为图间距）
        avail_w = max(1.0, cell_w - 2 * margin)
        avail_h = max(1.0, cell_h - 2 * margin)

        # 等比缩放以放入小格
        scale = min(avail_w / sw, avail_h / sh)

        new_w = sw * scale
        new_h = sh * scale

        # 将小页放到（r,c）格；坐标原点在左下
        x_left   = c * cell_w + (cell_w - new_w) / 2
        y_bottom = page_height - (r + 1) * cell_h + (cell_h - new_h) / 2

        tfm = (Transformation()
               .scale(sx=scale, sy=scale)
               .translate(tx=x_left, ty=y_bottom))
        big.merge_transformed_page(src, tfm)

    # 写出拼好的单页 PDF
    output_pdf = str(output_pdf)
    Path(output_pdf).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdf, "wb") as f:
        writer.write(f)

    # 同步导出 PNG（把拼好的一页渲染为位图）
    if output_png:
        doc = fitz.open(output_pdf)
        page = doc[0]
        mat = fitz.Matrix(render_zoom, render_zoom)
        pix = page.get_pixmap(matrix=mat, alpha=False)
        pix.save(output_png)
        doc.close()

    return output_pdf, output_png

def collage_png_pages_to_single(
    input_png: list[str],
    output_png: str,
    page_width: float = 1920,      # 16:9 画布宽（像素）
    page_height: float = 1080,     # 16:9 画布高
    cols: int | None = None,       # 固定列数；不设则自动按 16:9 估计
    margin: float = 2.0            # 每个小格四周留白；设更小=图间距更小
):
    from PIL import Image

    n = len(input_png)
    if n == 0:
        raise ValueError("输入 PNG 列表为空。")

    # 自动估列数，尽量匹配 16:9
    aspect = page_width / page_height
    if cols is None:
        cols = math.ceil(math.sqrt(n * aspect))
    rows = math.ceil(n / cols)

    cell_w = page_width / cols
    cell_h = page_height / rows

    big_img = Image.new("RGB", (int(page_width), int(page_height)), (255, 255, 255))

    for i, png_path in enumerate(input_png):
        r, c = divmod(i, cols)

        img = Image.open(png_path)
        sw, sh = img.size

        # 可用区域（留出 margin 作为图间距）
        avail_w = max(1.0, cell_w - 2 * margin)
        avail_h = max(1.0, cell_h - 2 * margin)

        # 等比缩放以放入小格
        scale = min(avail_w / sw, avail_h / sh)

        new_w = int(sw * scale)
        new_h = int(sh * scale)
        img_resized = img.resize((new_w, new_h), Image.Resampling.LANCZOS)

        # 将小图放到（r,c）格；坐标原点在左上
        x_left   = int(c * cell_w + (cell_w - new_w) / 2)
        y_top    = int(r * cell_h + (cell_h - new_h) / 2)

        big_img.paste(img_resized, (x_left, y_top))

    # 写出拼好的单张 PNG
    output_png = str(output_png)
    Path(output_png).parent.mkdir(parents=True, exist_ok=True)
    big_img.save(output_png)
    return output_png

input_path = "/home/wuct/MetaData/DATA/PbPb/apass4/ptShift/cutvar_woRef_combined/ptcenters_hist"
cuts = [cut for cut in Path(input_path).iterdir() if cut.is_dir() and 'ptcenter' in cut.name]
cuts.sort()
import re
pt = re.compile(r'pt_(\d+)_(\d+)')
for iCut, cut in enumerate(cuts):
    inpdf = []
    for file in cut.rglob("*.png"):
        if 'fPt' in file.name:
            inpdf.append(str(file))
    inpdf.sort(key=lambda x: int(pt.search(x).group(0).split('_')[1]))
    print(inpdf)
    outpdf = input_path + "/collaged_fPt" + f"_ptcenter{iCut}.pdf"
    outpng = outpdf.replace('.pdf', '.png')

    collage_png_pages_to_single(
        input_png=inpdf,
        output_png=outpng,
        page_width=1920, page_height=1080,  # 16:9
        cols=5,           # 让函数自动估列数；也可手动设 cols=4 之类
        margin=0,          # 图与图之间更紧凑；想更紧就设成 0
    )

# inpdf = "/home/wuct/ALICE/local/RTools/cfAnRes/tools/template/fitSpectrum_D0_PbPb_FT0C.pdf"
# outpdf = inpdf.replace('.pdf', '_collaged.pdf')
# outpng = inpdf.replace('.pdf', '_collaged.png')

# collage_pdf_pages_to_single(
#     input_pdf=inpdf,
#     output_pdf=outpdf,
#     output_png=outpng,   # 同时输出 PNG
#     page_width=1920, page_height=860,  # 16:9
#     cols=4,           # 让函数自动估列数；也可手动设 cols=4 之类
#     margin=0.5,          # 图与图之间更紧凑；想更紧就设成 0
#     render_zoom=3     # PNG 更清晰（文件更大），2~3 之间按需调
# )

