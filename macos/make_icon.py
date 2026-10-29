"""Draw the app icon (a codon wheel on a blue tile, one codon lit orange) and write AppIcon.icns."""
import math
import os
import subprocess
import sys
import tempfile

from PIL import Image, ImageDraw

CLASS_COLS = [(142, 107, 224), (31, 158, 137), (214, 64, 159), (47, 111, 217)]


def draw(size: int) -> Image.Image:
    s = size * 4  # supersample
    img = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    m = int(s * 0.06)
    tile = Image.new("RGBA", (s, s))
    td = ImageDraw.Draw(tile)
    for y in range(s):  # vertical gradient
        t = y / s
        c = (int(10 + 84 * t), int(132 - 40 * t), int(255 - 25 * t), 255)
        td.line((0, y, s, y), fill=c)
    mask = Image.new("L", (s, s), 0)
    ImageDraw.Draw(mask).rounded_rectangle((m, m, s - m, s - m), radius=int(s * 0.22), fill=255)
    img.paste(tile, (0, 0), mask)
    d = ImageDraw.Draw(img)
    cx = cy = s / 2
    r_out, r_mid, r_in, r_core = s * 0.36, s * 0.27, s * 0.19, s * 0.08
    # outer ring: 64 codon segments coloured by amino-acid class
    for i in range(64):
        a0, a1 = i * 360 / 64 - 90, (i + 1) * 360 / 64 - 90
        col = CLASS_COLS[(i * 7 + i // 5) % 4] + (235,)
        if i == 30:
            col = (255, 159, 10, 255)
        d.pieslice((cx - r_out, cy - r_out, cx + r_out, cy + r_out), a0 + 0.6, a1 - 0.6, fill=col)
    # the reassigned codon stands out from the ring
    r_hi = s * 0.405
    d.pieslice((cx - r_hi, cy - r_hi, cx + r_hi, cy + r_hi), 30 * 360 / 64 - 90 + 0.6, 31 * 360 / 64 - 90 - 0.6,
               fill=(255, 159, 10, 255))
    d.ellipse((cx - r_mid, cy - r_mid, cx + r_mid, cy + r_mid), fill=(255, 255, 255, 255))
    for i in range(16):
        a = math.radians(i * 360 / 16 - 90)
        d.line((cx + r_in * math.cos(a), cy + r_in * math.sin(a), cx + r_mid * math.cos(a), cy + r_mid * math.sin(a)),
               fill=(200, 205, 215, 255), width=max(1, s // 220))
    d.ellipse((cx - r_in, cy - r_in, cx + r_in, cy + r_in), fill=(245, 245, 247, 255), outline=(200, 205, 215, 255), width=max(1, s // 220))
    for i in range(4):
        a = math.radians(i * 90 - 90)
        d.line((cx + r_core * math.cos(a), cy + r_core * math.sin(a), cx + r_in * math.cos(a), cy + r_in * math.sin(a)),
               fill=(200, 205, 215, 255), width=max(1, s // 200))
    d.ellipse((cx - r_core, cy - r_core, cx + r_core, cy + r_core), fill=(10, 132, 255, 255))
    return img.resize((size, size), Image.LANCZOS)


def main(out: str) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        iconset = os.path.join(tmp, "AppIcon.iconset")
        os.mkdir(iconset)
        for base in (16, 32, 128, 256, 512):
            draw(base).save(os.path.join(iconset, f"icon_{base}x{base}.png"))
            draw(base * 2).save(os.path.join(iconset, f"icon_{base}x{base}@2x.png"))
        subprocess.run(["iconutil", "-c", "icns", iconset, "-o", out], check=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "AppIcon.icns")
