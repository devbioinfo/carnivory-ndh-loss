"""Original, schematic plant silhouettes (drawn here from simple shapes; no external artwork).
Each icon is drawn in a unit box [0,1]x[0,1] inside its own inset axes."""
import numpy as np
from matplotlib.patches import Ellipse, Polygon, Circle, PathPatch, FancyBboxPatch, Wedge
from matplotlib.path import Path
def _ax(fig, rect):
    a = fig.add_axes(rect); a.set_xlim(0, 1); a.set_ylim(0, 1); a.axis('off'); a.set_aspect('equal', adjustable='datalim'); return a
def _stem(a, pts, c, lw=1.0): a.plot(*np.array(pts).T, color=c, lw=lw, solid_capstyle='round')
def _leaf(a, x, y, L, W, ang, c, z=2):
    a.add_patch(Ellipse((x + 0.5 * L * np.cos(np.radians(ang)), y + 0.5 * L * np.sin(np.radians(ang))), L, W, angle=ang, fc=c, ec='none', zorder=z))
def ground(a, c='#c3c2b7'): a.plot([0.05, 0.95], [0.04, 0.04], color=c, lw=0.6)
def autotroph(a, c):
    _stem(a, [(0.5, 0.04), (0.5, 0.78)], c, 1.1)
    for y, s in [(0.25, 1), (0.42, -1), (0.58, 1), (0.72, -1)]:
        _leaf(a, 0.5, y, 0.32, 0.12, 90 - s * 58, c)
    _leaf(a, 0.5, 0.78, 0.2, 0.1, 90, c); ground(a)
def sundew(a, c):  # rosette of spoon leaves with glandular tentacles (Drosera-like)
    for ang in np.linspace(15, 165, 6):
        r = np.radians(ang); x1, y1 = 0.5 + 0.28 * np.cos(r), 0.06 + 0.42 * np.sin(r)
        _stem(a, [(0.5, 0.06), (x1, y1)], c, 0.8)
        a.add_patch(Circle((x1 + 0.06 * np.cos(r), y1 + 0.06 * np.sin(r)), 0.07, fc=c, ec='none'))
        for k in np.linspace(-1, 1, 5):
            rr = r + k * 0.9; xe, ye = x1 + 0.06 * np.cos(r) + 0.11 * np.cos(rr), y1 + 0.06 * np.sin(r) + 0.11 * np.sin(rr)
            a.plot([x1 + 0.06 * np.cos(r), xe], [y1 + 0.06 * np.sin(r), ye], color=c, lw=0.35); a.add_patch(Circle((xe, ye), 0.018, fc=c, ec='none'))
    _stem(a, [(0.5, 0.06), (0.52, 0.95)], c, 0.7); a.add_patch(Circle((0.53, 0.95), 0.03, fc=c, ec='none')); ground(a)
def trumpet(a, c):  # tall tubular pitchers with hood (Sarracenia-like)
    for x0, h, s in [(0.32, 0.86, 0.9), (0.55, 0.95, 1.0), (0.74, 0.7, 0.8)]:
        w = 0.07 * s; poly = [(x0 - 0.02, 0.04), (x0 - w, h), (x0 + w, h), (x0 + 0.02, 0.04)]
        a.add_patch(Polygon(poly, fc=c, ec='none'))
        a.add_patch(Ellipse((x0 + 0.02, h + 0.02), 2.4 * w, 0.09, angle=-15, fc=c, ec='none'))
    ground(a)
def nepenthes(a, c):  # hanging pitcher on a tendril with lid (Nepenthes-like)
    _stem(a, [(0.08, 0.98), (0.2, 0.86)], c, 1.0)
    _leaf(a, 0.2, 0.86, 0.48, 0.13, -12, c)
    a.plot([0.66, 0.76, 0.74, 0.66], [0.76, 0.7, 0.58, 0.53], color=c, lw=0.8)
    t = np.linspace(0, 2 * np.pi, 80); xs = 0.6 + 0.13 * np.sin(t) * (1 + 0.35 * np.cos(t)); ys = 0.3 - 0.24 * np.cos(t) * 0.9
    a.add_patch(Polygon(np.c_[xs, ys], fc=c, ec='none'))
    a.add_patch(Ellipse((0.6, 0.5), 0.2, 0.05, fc='white', ec=c, lw=0.6)); a.add_patch(Ellipse((0.68, 0.56), 0.18, 0.06, angle=25, fc=c, ec='none'))
def cephalotus(a, c):  # squat lidded pitchers with ribs
    for x0, s in [(0.3, 0.8), (0.65, 1.0)]:
        a.add_patch(Ellipse((x0, 0.3 * s), 0.3 * s, 0.42 * s, fc=c, ec='none'))
        for dx in (-0.06, 0.06): a.plot([x0 + dx * s, x0 + dx * s * 0.6], [0.14 * s, 0.48 * s], color='white', lw=0.5)
        a.add_patch(Ellipse((x0 + 0.04 * s, 0.56 * s), 0.3 * s, 0.1 * s, angle=20, fc=c, ec='none'))
    ground(a)
def bladderwort(a, c):  # emergent flower stalk, submerged stolon with bladders (Utricularia-like)
    a.axhspan(0, 0.42, color='#e6eef8', lw=0); _stem(a, [(0.5, 0.42), (0.5, 0.9)], c, 0.9)
    a.add_patch(Ellipse((0.53, 0.92), 0.16, 0.1, fc=c, ec='none')); a.add_patch(Ellipse((0.5, 0.83), 0.1, 0.08, fc=c, ec='none'))
    _stem(a, [(0.08, 0.3), (0.92, 0.3)], c, 0.8)
    for x in np.linspace(0.14, 0.86, 6):
        _stem(a, [(x, 0.3), (x + 0.04, 0.18)], c, 0.5); a.add_patch(Circle((x + 0.05, 0.15), 0.035, fc=c, ec='none'))
        _stem(a, [(x, 0.3), (x - 0.05, 0.38)], c, 0.4)
def byblis(a, c):  # erect, very narrow glandular leaves
    for ang in (60, 75, 90, 105, 120):
        r = np.radians(ang); x1, y1 = 0.5 + 0.45 * np.cos(r), 0.04 + 0.85 * np.sin(r)
        _stem(a, [(0.5, 0.04), (x1, y1)], c, 0.9)
        for f in np.linspace(0.25, 0.95, 6):
            xp, yp = 0.5 + (x1 - 0.5) * f, 0.04 + (y1 - 0.04) * f
            a.add_patch(Circle((xp + 0.025, yp), 0.012, fc=c, ec='none')); a.add_patch(Circle((xp - 0.025, yp), 0.012, fc=c, ec='none'))
    ground(a)
def roridula(a, c):  # small shrub, branched, sticky linear leaves at tips
    _stem(a, [(0.5, 0.04), (0.5, 0.5)], c, 1.2)
    for x2, y2 in [(0.25, 0.75), (0.5, 0.85), (0.75, 0.72)]:
        _stem(a, [(0.5, 0.45), (x2, y2)], c, 0.9)
        for ang in np.linspace(0, 330, 10):
            r = np.radians(ang); a.plot([x2, x2 + 0.11 * np.cos(r)], [y2, y2 + 0.11 * np.sin(r)], color=c, lw=0.5)
    ground(a)
def triantha(a, c):  # grass-like leaves, sticky glandular raceme
    for ang in (60, 78, 100, 120): _leaf(a, 0.5, 0.04, 0.42, 0.08, ang, c)
    _stem(a, [(0.52, 0.04), (0.54, 0.95)], c, 0.9)
    for y in np.linspace(0.5, 0.92, 6):
        for dx in (-0.06, 0.06): a.add_patch(Circle((0.54 + dx, y), 0.045, fc=c, ec='none'))
    ground(a)
def tank(a, c):  # tubular tank rosette (Brocchinia-like)
    for ang, L in [(70, 0.85), (80, 0.95), (90, 1.0), (100, 0.95), (110, 0.85)]:
        _leaf(a, 0.5, 0.04, 0.8 * L, 0.14, ang, c)
    ground(a)
def broomrape(a, c):  # leafless spike with scale leaves on a host root (Orobanche-like)
    a.plot([0.05, 0.95], [0.08, 0.1], color='#c3c2b7', lw=1.6)  # host root
    a.add_patch(Ellipse((0.5, 0.11), 0.16, 0.08, fc=c, ec='none'))  # haustorial tuber
    a.add_patch(Polygon([(0.44, 0.12), (0.46, 0.92), (0.54, 0.92), (0.56, 0.12)], fc=c, ec='none'))
    for y in np.linspace(0.4, 0.88, 6):
        for s in (-1, 1): a.add_patch(Ellipse((0.5 + s * 0.08, y), 0.12, 0.06, angle=s * 35, fc=c, ec='none'))
    for y in np.linspace(0.18, 0.36, 3):
        for s in (-1, 1): a.add_patch(Ellipse((0.5 + s * 0.05, y), 0.07, 0.025, angle=s * 50, fc=c, ec='none'))
def hemiparasite(a, c):  # green leafy plant attached to host root by haustoria (Pedicularis/Striga-like)
    a.plot([0.05, 0.95], [0.08, 0.1], color='#c3c2b7', lw=1.6)
    _stem(a, [(0.5, 0.1), (0.5, 0.88)], c, 1.0)
    for y, s in [(0.3, 1), (0.42, -1), (0.55, 1), (0.66, -1)]: _leaf(a, 0.5, y, 0.26, 0.07, 90 - s * 62, c)
    for y in (0.76, 0.84, 0.92): a.add_patch(Circle((0.5, y), 0.04, fc=c, ec='none'))
    _stem(a, [(0.5, 0.1), (0.35, 0.09)], c, 0.6); _stem(a, [(0.5, 0.1), (0.66, 0.1)], c, 0.6)
def monotropa(a, c):  # pale leafless stems hooked at the tip, nodding bell flower (Monotropa-like)
    for x0, h in [(0.32, 0.72), (0.6, 0.85)]:
        t = np.linspace(0, 1, 40); xs = x0 + 0.14 * (1 - np.cos(t * np.pi / 2)); ys = 0.04 + h * np.sin(t * np.pi * 0.62) / np.sin(np.pi * 0.62 * 0.81)
        ys = 0.04 + h * (np.sin(t * np.pi * 0.62))
        a.plot(xs, ys, color=c, lw=2.0, solid_capstyle='round')
        xe, ye = xs[-1], ys[-1]
        a.add_patch(Polygon([(xe - 0.035, ye), (xe - 0.07, ye - 0.2), (xe + 0.07, ye - 0.2), (xe + 0.035, ye)], fc=c, ec='none'))
        for f in (0.2, 0.4, 0.6): a.add_patch(Ellipse((x0 + 0.14 * (1 - np.cos(f * np.pi / 2)) + 0.03, 0.04 + h * np.sin(f * np.pi * 0.62)), 0.05, 0.02, angle=40, fc=c, ec='none'))
    ground(a)
KIND = dict(autotroph=autotroph, sundew=sundew, trumpet=trumpet, nepenthes=nepenthes, cephalotus=cephalotus, bladderwort=bladderwort,
            byblis=byblis, roridula=roridula, triantha=triantha, tank=tank, broomrape=broomrape, hemiparasite=hemiparasite, monotropa=monotropa)
def icon(fig, kind, rect, color):
    a = _ax(fig, rect); KIND[kind](a, color); return a
def icon_in(ax, kind, bounds, color):
    a = ax.inset_axes(bounds); a.set_xlim(0, 1); a.set_ylim(0, 1); a.axis('off'); a.set_aspect('equal', adjustable='datalim'); KIND[kind](a, color); return a
def succulent(a, c):  # ribbed columnar succulent with areoles (CAM)
    a.add_patch(FancyBboxPatch((0.36, 0.04), 0.28, 0.8, boxstyle='round,pad=0,rounding_size=0.14', fc=c, ec='none'))
    a.add_patch(FancyBboxPatch((0.12, 0.3), 0.16, 0.36, boxstyle='round,pad=0,rounding_size=0.08', fc=c, ec='none')); _stem(a, [(0.2, 0.32), (0.2, 0.26), (0.38, 0.26)], c, 2.5)
    for x in (0.44, 0.5, 0.56): a.plot([x, x], [0.1, 0.78], color='white', lw=0.4)
    for y in np.linspace(0.15, 0.75, 6):
        for x in (0.36, 0.64): a.plot([x, x + (0.05 if x > 0.5 else -0.05)], [y, y + 0.02], color=c, lw=0.4)
    ground(a)
KIND['succulent'] = succulent
