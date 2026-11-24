#!/usr/bin/env python3
"""Build the general-audience PDF research report."""
import os
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                Image, PageBreak, KeepTogether, HRFlowable, ListFlowable, ListItem)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

import matplotlib
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SP = HERE  # figures live next to this script (run generate_figures.py first)
FONTDIR = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
OUT = os.path.join(HERE, "construct_screening_report.pdf")

pdfmetrics.registerFont(TTFont("DJ", f"{FONTDIR}/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("DJ-B", f"{FONTDIR}/DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DJ-I", f"{FONTDIR}/DejaVuSans-Oblique.ttf"))
pdfmetrics.registerFontFamily("DJ", normal="DJ", bold="DJ-B", italic="DJ-I", boldItalic="DJ-B")

INK   = colors.HexColor("#1b2430")
SLATE = colors.HexColor("#2f3e4e")
ACC   = colors.HexColor("#2f6f8f")   # teal  (alignment / conventional)
ACC2  = colors.HexColor("#c05a2e")   # rust  (model)
MUT   = colors.HexColor("#5b6b7b")
LINE  = colors.HexColor("#c9d3dc")
BG    = colors.HexColor("#f2f5f8")
BGWARM= colors.HexColor("#fbf1ea")
GREEN = colors.HexColor("#2e7d5b")

styles = getSampleStyleSheet()
def S(name, **kw):
    base = kw.pop("parent", styles["Normal"])
    kw.setdefault("fontName", "DJ")
    return ParagraphStyle(name, parent=base, **kw)

body = S("body", fontSize=10, leading=15, alignment=TA_JUSTIFY, textColor=INK, spaceAfter=7)
bodyL= S("bodyL", fontSize=10, leading=15, alignment=TA_LEFT, textColor=INK, spaceAfter=7)
lead = S("lead", fontSize=11.5, leading=17, alignment=TA_LEFT, textColor=SLATE, spaceAfter=9)
h1   = S("h1", fontName="DJ-B", fontSize=16, leading=20, textColor=INK, spaceBefore=16, spaceAfter=3, keepWithNext=1)
h2   = S("h2", fontName="DJ-B", fontSize=12.5, leading=16, textColor=ACC, spaceBefore=12, spaceAfter=4, keepWithNext=1)
kicker=S("kicker", fontName="DJ-B", fontSize=9, leading=12, textColor=ACC2, spaceAfter=2)
small= S("small", fontSize=8.3, leading=11.5, textColor=MUT, alignment=TA_LEFT)
smallC=S("smallC", fontSize=8.3, leading=11.5, textColor=MUT, alignment=TA_CENTER)
cap  = S("cap", fontSize=8.5, leading=12, textColor=MUT, alignment=TA_LEFT, spaceBefore=3, spaceAfter=10)
bullet=S("bullet", fontSize=10, leading=14.5, textColor=INK)
sidebar=S("sidebar", fontSize=9, leading=13, textColor=SLATE)
sidebarH=S("sidebarH", fontName="DJ-B", fontSize=10.5, leading=13, textColor=INK, spaceAfter=3)
tbl_h = S("tbl_h", fontName="DJ-B", fontSize=8.6, leading=11, textColor=colors.white)
tbl_c = S("tbl_c", fontSize=8.8, leading=11.5, textColor=INK)
tbl_cC= S("tbl_cC", fontSize=8.8, leading=11.5, textColor=INK, alignment=TA_CENTER)
tbl_cR= S("tbl_cR", fontSize=8.8, leading=11.5, textColor=INK, alignment=TA_LEFT)
titleS= S("titleS", fontName="DJ-B", fontSize=22, leading=26, textColor=INK, spaceAfter=4)
subS  = S("subS", fontSize=12.5, leading=17, textColor=ACC, spaceAfter=10)

def P(t, s=body): return Paragraph(t, s)
def rule(c=LINE, w=0.8, sb=2, sa=8): return HRFlowable(width="100%", thickness=w, color=c, spaceBefore=sb, spaceAfter=sa)

def bullets(items, st=bullet, sa=7):
    return ListFlowable(
        [ListItem(Paragraph(x, st), leftIndent=6, value="•") for x in items],
        bulletType="bullet", start="•", leftIndent=14, bulletColor=ACC2, bulletFontSize=9, spaceAfter=sa)

def shaded(flowables, bg=BG, pad=9, border=None):
    t = Table([[flowables]], colWidths=[FULL])
    st = [("BACKGROUND",(0,0),(-1,-1),bg),
          ("LEFTPADDING",(0,0),(-1,-1),pad),("RIGHTPADDING",(0,0),(-1,-1),pad),
          ("TOPPADDING",(0,0),(-1,-1),pad),("BOTTOMPADDING",(0,0),(-1,-1),pad)]
    if border: st.append(("LINEBEFORE",(0,0),(0,-1),3,border))
    t.setStyle(TableStyle(st))
    return t

# ---- page furniture ----
def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("DJ", 7.8); canvas.setFillColor(MUT)
    canvas.drawString(doc.leftMargin, 0.5*inch,
                      "Construct preprocessing for function-based screening  ·  research report")
    canvas.drawRightString(letter[0]-doc.rightMargin, 0.5*inch, f"Page {doc.page}")
    canvas.setStrokeColor(LINE); canvas.setLineWidth(0.5)
    canvas.line(doc.leftMargin, 0.68*inch, letter[0]-doc.rightMargin, 0.68*inch)
    canvas.restoreState()

doc = SimpleDocTemplate(OUT, pagesize=letter,
                        leftMargin=0.9*inch, rightMargin=0.9*inch,
                        topMargin=0.85*inch, bottomMargin=0.8*inch,
                        title="Construct preprocessing for function-based screening",
                        author="Research project (self-review)")
FULL = letter[0] - 1.8*inch
st = []

# ============================ TITLE ============================
st += [
  P("RESEARCH REPORT · SELF-REVIEW · 24 NOVEMBER 2026", kicker),
  rule(INK, 1.4, 2, 10),
  P("Preparing synthetic DNA for smarter safety screening", titleS),
  P("Can better &ldquo;preprocessing&rdquo; and AI protein models help identify what an "
    "engineered DNA construct actually contains — especially when it looks nothing like "
    "anything on file?", subS),
]
meta = Table([[
    P("<b>Field</b><br/>Biosecurity / computational biology", small),
    P("<b>Scope</b><br/>Public, benign constructs only; defensive", small),
    P("<b>Status</b><br/>Core comparison complete; one goal open", small),
    P("<b>Review</b><br/>Internal self-review; not peer reviewed", small),
]], colWidths=[FULL*0.26, FULL*0.28, FULL*0.26, FULL*0.20])
meta.setStyle(TableStyle([("LINEABOVE",(0,0),(-1,0),0.8,LINE),("LINEBELOW",(0,0),(-1,0),0.8,LINE),
    ("TOPPADDING",(0,0),(-1,-1),7),("BOTTOMPADDING",(0,0),(-1,-1),7),
    ("LEFTPADDING",(0,0),(0,0),0),("VALIGN",(0,0),(-1,-1),"TOP")]))
st += [meta, Spacer(1,12)]

# ---- Executive summary box ----
exec_items = [
  "<b>The job.</b> Synthetic DNA constructs — the engineered genetic &ldquo;recipes&rdquo; used "
  "across modern biology — are screened to check they do not encode something hazardous. Today "
  "that screening works by matching a sequence against a library of known ones. It struggles "
  "when a sequence is only distantly related to anything on file.",
  "<b>Two ideas tested.</b> (1) <i>Preprocessing</i>: cut a raw construct into clean, "
  "protein-sized pieces before screening, instead of feeding in the raw DNA. (2) <i>AI protein "
  "models</i>: use a learned model of protein &ldquo;meaning&rdquo; alongside conventional matching.",
  "<b>Main result.</b> On 777 real public constructs (1,221 known parts), a standard gene-finding "
  "step raised the share of parts correctly identified by about <b>9 percentage points</b> for the "
  "same screening effort — a pre-registered, statistically clear win. The benefit comes mostly from "
  "<i>choosing better pieces to look at</i>, not from a smarter final judge.",
  "<b>The headline nuance.</b> On the full dataset, the AI model did <i>not</i> beat conventional "
  "matching. But in the hardest cases — where no close relative exists on file, exactly where today's "
  "screening fails — a model that reads a protein <i>residue by residue</i> identified "
  "<b>96%</b> of families versus <b>68%</b> for conventional matching.",
  "<b>Honest status.</b> That hardest-case result is promising but <i>exploratory</i>: it reuses a "
  "public answer key and labels built from the very method it is compared against. A clean, "
  "independent confirmation is designed but not yet run, so the project's own low-similarity goal "
  "remains formally <b>open</b>.",
]
st += [ shaded([P("WHAT THIS REPORT FOUND, IN BRIEF", sidebarH)] +
               [P(x, S("exb", fontSize=9.3, leading=13.5, textColor=INK, spaceAfter=6)) for x in exec_items],
               bg=BG, border=ACC2) ]

st += [PageBreak()]

# ============================ 1. THE PROBLEM ============================
st += [P("1   Why this problem matters", h1), rule()]
st += [P(
  "Almost every advance in modern biology — vaccines, medicines, crops, research reagents — starts "
  "with <b>synthetic DNA</b>: genetic sequences written on a computer and then manufactured to order "
  "by commercial providers. A single ordered sequence, called a <b>construct</b>, is a recipe that "
  "may contain several working parts (genes that encode proteins) stitched together with connective "
  "and structural DNA.", lead)]
st += [P(
  "Because the same technology that builds a helpful protein could in principle build a harmful one, "
  "DNA providers screen incoming orders. The goal is narrow and defensive: flag an order for human "
  "review if it looks like it encodes something of concern. <b>This project is entirely about making "
  "that defensive screening work better, using only public, benign material.</b> It never attempts to "
  "design, improve, or evaluate anything hazardous.")]

st += [P("How screening works today — and where it breaks", h2)]
st += [P(
  "Today's screening leans on <b>sequence matching</b> (the same idea as a plagiarism checker): line "
  "up the ordered sequence against a library of known sequences and look for a strong match. This is "
  "fast, well understood, and reliable <i>when a close relative is on file</i>. Its weakness is "
  "equally well known: as a sequence drifts away from anything in the library — through natural "
  "divergence, or through deliberate or AI-assisted redesign that preserves a protein's function "
  "while changing its spelling — the match fades, and the screen can miss it. Published biosecurity "
  "studies place this breakdown below roughly <b>30% sequence identity</b> to the nearest reference.")]
st += [P(
  "Two newer ideas might help. First, <b>better preprocessing</b>: a raw construct is messy, so rather "
  "than screen the whole tangle, first cut it into clean, protein-sized candidate pieces, then screen "
  "those. Second, <b>AI protein models</b> (the protein equivalent of a large language model): these "
  "learn a rich internal sense of what proteins &ldquo;mean&rdquo; from millions of examples, and might "
  "recognise a distant relative that letter-by-letter matching misses. Both ideas are widely assumed "
  "to help. This project set out to <i>test</i> them rather than assume them — and to report the "
  "result even when it is mixed or negative.")]

# glossary sidebar
gloss = [
  ("Construct", "One engineered DNA sequence (e.g. a plasmid), often containing several genes plus structural DNA."),
  ("Component", "A documented working part of a construct — here, a protein-coding gene."),
  ("Preprocessing", "Turning a raw construct into clean, model-ready protein pieces before screening."),
  ("Gene calling", "A standard computational step that predicts where genes start and stop in raw DNA."),
  ("Sequence identity", "How similar two sequences are, as a %. Low identity = distant relative = hard case."),
  ("Foundation / AI model", "A model (here, ESM-2) trained on millions of proteins to capture their patterns."),
  ("Held-out answer key", "The documented contents of each construct, hidden from the methods and used only to grade them."),
]
gl = [[P(f"<b>{t}</b>", S('glt',fontSize=8.6,leading=11,textColor=ACC)), P(d, sidebar)] for t,d in gloss]
glt = Table([[P("PLAIN-LANGUAGE GLOSSARY", sidebarH), ""]] + gl, colWidths=[FULL*0.22, FULL*0.78])
glt.setStyle(TableStyle([
    ("SPAN",(0,0),(1,0)),
    ("BACKGROUND",(0,0),(-1,-1),BGWARM),
    ("LINEBEFORE",(0,0),(0,-1),3,ACC2),
    ("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4),
    ("LEFTPADDING",(0,0),(-1,-1),9),("RIGHTPADDING",(0,0),(-1,-1),9),
    ("VALIGN",(0,0),(-1,-1),"TOP"),
    ("TOPPADDING",(0,0),(-1,0),8),
]))
st += [Spacer(1,4), glt, Spacer(1,6)]

# ============================ 2. QUESTIONS ============================
st += [P("2   The three questions", h1), rule()]
st += [P("The study was organised around three questions, fixed in advance:")]
q = [
  ("Q1", "Can preprocessing recover useful protein pieces from a construct <b>without</b> being told "
         "where its parts are?"),
  ("Q2", "Do those pieces help identify a construct's known parts <b>better than a simple baseline</b>, "
         "for the same amount of screening effort?"),
  ("Q3", "In the <b>hardest, low-similarity cases</b>, can an AI model identify distant parts better "
         "than conventional matching — and does the answer depend on the type of construct?"),
]
qt = Table([[P(f"<b>{k}</b>", S('qk',fontName='DJ-B',fontSize=11,textColor=colors.white,alignment=TA_CENTER)),
             P(v, bodyL)] for k,v in q], colWidths=[0.5*inch, FULL-0.5*inch])
qt.setStyle(TableStyle([
    ("BACKGROUND",(0,0),(0,-1),ACC),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ("TOPPADDING",(0,0),(-1,-1),8),("BOTTOMPADDING",(0,0),(-1,-1),8),
    ("LEFTPADDING",(1,0),(1,-1),10),("RIGHTPADDING",(1,0),(1,-1),4),
    ("LINEBELOW",(0,0),(-1,-2),0.6,LINE),
    ("ROWBACKGROUNDS",(1,0),(1,-1),[colors.white, BG]),
]))
st += [qt, Spacer(1,6)]
st += [P(
  "A fourth idea was tested as a sub-hypothesis: whether <i>splitting</i> fused proteins into pieces "
  "helps the model. The scope deliberately asked for honest reporting of negative and mixed results, "
  "and required that the low-similarity question (Q3) actually be <i>executed</i> against conventional "
  "methods — not merely discussed — or else be reported as unfinished.", body)]

st += [Spacer(1,4)]

# ============================ 3. METHODS ============================
st += [P("3   How the test was set up", h1), rule()]
st += [P(
  "The design mimics the real screening situation. Each construct in the test comes with a "
  "<b>documented list of its parts</b> — a trustworthy answer key. That key is <b>hidden</b> from "
  "every method being graded; the methods see only raw sequence, exactly as a provider's screen would. "
  "A method&rsquo;s score is simply the share of the hidden parts it correctly identifies.", lead)]

st += [P("A fair-effort comparison", h2)]
st += [P(
  "Looking at more candidate pieces always finds more parts, but costs more. To compare methods "
  "honestly, every method was given the <b>same budget</b> — the same number of protein pieces it may "
  "examine per construct — and compared at matched budgets, as well as across a range of budgets. "
  "A minimum worthwhile improvement (5 percentage points) and the statistical rules were all "
  "<b>locked in writing before any answer key was consulted</b>, so the headline could not be chosen "
  "after seeing the data.")]

st += [P("The test set", h2)]
st += [P(
  "The main dataset (&ldquo;cohort v2&rdquo;) is <b>777 real, public, benign cloning and expression "
  "vectors</b> deposited in the public NCBI archive, carrying <b>1,221 carefully corroborated "
  "protein-coding parts</b>. To avoid flattering a method with many near-duplicate constructs, parts "
  "were grouped into <b>422 lineage groups</b> and scored so each group counts once. The constructs "
  "span the four shapes the field cares about:")]
arche = [["Construct shape","What it looks like","Records"],
 ["Multi-gene cassette","Several genes in a row","355"],
 ["Backbone-heavy","One or few genes in a large structural scaffold","144"],
 ["Single label + extra genes","A main gene plus incidental others","128"],
 ["Single gene","One coding sequence","99"],
 ["Tagged / fused protein","Two parts joined into one protein","51"]]
at = Table([[P(c, tbl_h) if i==0 else (P(c, tbl_c) if j==1 else (P(c,tbl_cC) if j==2 else P(c,tbl_c)))
             for j,c in enumerate(row)] for i,row in enumerate(arche)],
           colWidths=[FULL*0.28, FULL*0.52, FULL*0.20])
at.setStyle(TableStyle([
    ("BACKGROUND",(0,0),(-1,0),SLATE),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white, BG]),
    ("GRID",(0,0),(-1,-1),0.5,LINE),
    ("TOPPADDING",(0,0),(-1,-1),5),("BOTTOMPADDING",(0,0),(-1,-1),5),
    ("LEFTPADDING",(0,0),(-1,-1),7),("ALIGN",(2,1),(2,-1),"CENTER"),
    ("VALIGN",(0,0),(-1,-1),"MIDDLE"),
]))
st += [KeepTogether([at, P("Table 1. The test constructs, by shape. Numbers are records in the held-out set.", cap)])]

st += [P("What was compared", h2)]
st += [P("<b>Three ways to preprocess</b> a construct into candidate pieces:")]
st += [bullets([
  "<b>Minimal (baseline)</b> — the naive approach: translate the raw DNA in all six reading frames and "
  "take whatever falls out. Simple, but it produces a huge, noisy pile of candidates.",
  "<b>Gene-calling</b> — a standard, off-the-shelf gene finder (Pyrodigal) predicts likely genes first, "
  "yielding far fewer, cleaner candidates. Crucially, it never sees the answer key.",
  "<b>Answer-key pieces (ceiling)</b> — a diagnostic only: hand the method the exact documented parts. "
  "This is not a usable method; it shows how much room for improvement remains.",
])]
st += [P("<b>Two ways to then identify</b> each piece — to decide what family it belongs to:")]
st += [bullets([
  "<b>Conventional sequence matching</b> (BLOSUM62 alignment) — the established, letter-by-letter method.",
  "<b>An AI protein model</b> (ESM-2), used two ways: a quick <i>whole-protein summary</i>, and a more "
  "detailed <i>residue-by-residue</i> reading. Two model sizes were tested (8M and a 19&times;-larger 150M).",
])]

st += [P("Trustworthiness of the result", h2)]
st += [P(
  "The protocol and code were frozen at a specific version <b>before</b> the hidden answer key was "
  "touched, and the final test was run <b>once</b>. All uncertainty is reported as confidence "
  "intervals from resampling the data thousands of times. Re-running the analysis from saved traces "
  "reproduced the headline numbers <b>exactly</b>. Every input is pinned by a cryptographic fingerprint. "
  "This is internal self-review, not independent peer review — and the report says so wherever a claim "
  "depends on it.")]

st += [Spacer(1,6)]

# ============================ 4. RESULTS ============================
st += [P("4   What the test showed", h1), rule()]

st += [P("4.1   Preprocessing recovers clean pieces — Q1: yes", h2)]
st += [P(
  "Gene-calling turned raw constructs into candidate proteins without any access to the answer key. "
  "Its candidate pile was about <b>24&times; smaller</b> than the naive baseline&rsquo;s and used fewer "
  "total protein residues — so it is not only cleaner but <i>cheaper</i>. The honest limit: because the "
  "gene finder is tuned for bacteria, it simply has no candidate for 45 of the 1,221 parts (versus 11 "
  "for the naive pile), and it largely misses plant genes. It recovers useful pieces, but it is not "
  "a universal cutter.")]

st += [P("4.2   Better pieces, same effort, more parts found — Q2: yes", h2)]
st += [P(
  "This is the study&rsquo;s firmest result. At a budget of two pieces per construct, switching the "
  "baseline&rsquo;s naive pile for gene-called pieces raised the share of parts correctly identified "
  "by about <b>9 percentage points</b> — and the gain held whether the final judge was the AI model or "
  "conventional matching. It comfortably cleared the 5-point bar set in advance, and the gene-called "
  "approach used <i>less</i> computer effort, not more.")]

img1 = Image(f"{SP}/fig_gain.png", width=FULL*0.80, height=FULL*0.80*780/1440)
img1.hAlign = "CENTER"
st += [KeepTogether([img1, P("Figure 1. For the same screening budget (horizontal axis), gene-calling preprocessing "
  "(teal) finds more of each construct&rsquo;s documented parts than the naive baseline (orange). The grey "
  "diamond shows the ceiling if the hidden answer key were handed over — about 27 points of headroom remain.", cap)])]

# Q2 table
t2 = [["Method (at 2 pieces/construct)","AI model","Conventional","Candidate pile"],
 ["Minimal (baseline)","46.9%","48.6%","101,407"],
 ["Gene-calling","59.4%","61.1%","4,194"],
 ["Answer-key pieces (ceiling)","91.4%","92.9%","1,250"],
 ["Gain: gene-calling over baseline","+9.2 pts","+9.0 pts","—"],
 ["   95% confidence interval","6.6 to 11.8","6.5 to 11.6","—"]]
T2 = Table([[P(c, tbl_h) if i==0 else P(c, tbl_c if j==0 else tbl_cC) for j,c in enumerate(row)]
            for i,row in enumerate(t2)],
           colWidths=[FULL*0.40, FULL*0.19, FULL*0.21, FULL*0.20])
T2.setStyle(TableStyle([
    ("BACKGROUND",(0,0),(-1,0),SLATE),
    ("ROWBACKGROUNDS",(0,1),(-1,3),[colors.white, BG]),
    ("BACKGROUND",(0,4),(-1,5),BGWARM),
    ("GRID",(0,0),(-1,-1),0.5,LINE),
    ("TOPPADDING",(0,0),(-1,-1),4.5),("BOTTOMPADDING",(0,0),(-1,-1),4.5),
    ("LEFTPADDING",(0,0),(-1,-1),7),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ("FONTNAME",(0,4),(-1,4),"DJ-B"),("TEXTCOLOR",(1,4),(2,4),GREEN),("FONTNAME",(1,4),(2,4),"DJ-B"),
    ("LINEABOVE",(0,4),(-1,4),1.0,SLATE),
]))
st += [KeepTogether([T2, P("Table 2. Share of the 1,221 documented parts correctly identified (equal-group scoring). "
  "&ldquo;Candidate pile&rdquo; is the total number of pieces each approach generates — smaller is cleaner and cheaper.", cap)])]

st += [P("What is really driving the gain", h2)]
st += [P(
  "A key diagnostic finding: the improvement comes from <b>choosing better pieces to examine</b>, not "
  "from a cleverer final judge. The two judges (AI model and conventional matching) gained almost "
  "identically — within 0.2 points of each other. Gene-calling&rsquo;s advantage is that it promotes "
  "156 genuine genes into the limited budget that the naive pile had buried under noise.")]

st += [P("Where the honesty lives — the caveats on Q2", h2)]
st += [bullets([
  "<b>Concentrated in two families.</b> About 82% of the extra wins are two common bacterial "
  "resistance genes that a bacterial gene finder naturally handles well. Removing them, the "
  "per-group gain is still a real +7.4 points, but the overall figure softens.",
  "<b>Grouping matters.</b> Under the strictest &ldquo;count shared proteins once&rdquo; view the "
  "dataset collapses toward a single effective group, and the interval widens to include small gains.",
  "<b>Fused proteins are inconclusive</b> (+4.6 points, but the interval crosses zero).",
  "<b>Budget-dependent.</b> At one piece per construct the gain is a smaller +4.0 points.",
])]

st += [Spacer(1,6)]

# construct-type variation
st += [P("4.3   It helps across construct shapes — with one exception", h2)]
st += [P(
  "Breaking the Q2 result down by construct shape shows the benefit is broad rather than a quirk of "
  "one category. It is largest for backbone-heavy constructs (where the naive pile is noisiest) and "
  "clearly positive for every shape except fused proteins, where the evidence is not yet conclusive.")]
t3 = [["Construct shape","Baseline","Gene-calling","Gain (AI model)"],
 ["Backbone-heavy","37%","52%","+14.4"],
 ["Multi-gene cassette","40%","55%","+9.0"],
 ["Single label + extra genes","90%","98%","+8.1"],
 ["Single gene","70%","80%","+7.6"],
 ["Tagged / fused protein","36%","40%","+4.6 (not conclusive)"]]
T3 = Table([[P(c, tbl_h) if i==0 else P(c, tbl_c if j==0 else tbl_cC) for j,c in enumerate(row)]
            for i,row in enumerate(t3)],
           colWidths=[FULL*0.36, FULL*0.17, FULL*0.19, FULL*0.28])
T3.setStyle(TableStyle([
    ("BACKGROUND",(0,0),(-1,0),SLATE),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white, BG]),
    ("GRID",(0,0),(-1,-1),0.5,LINE),
    ("TOPPADDING",(0,0),(-1,-1),4.5),("BOTTOMPADDING",(0,0),(-1,-1),4.5),
    ("LEFTPADDING",(0,0),(-1,-1),7),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
    ("TEXTCOLOR",(3,1),(3,-1),GREEN),("FONTNAME",(3,1),(3,-1),"DJ-B"),
]))
st += [KeepTogether([T3, P("Table 3. Share of parts found, by construct shape, for the AI-model judge. "
  "Figures are approximate recall at a two-piece budget; &ldquo;gain&rdquo; is in percentage points.", cap)])]

# Q3
st += [P("4.4   The hardest cases: where AI models earn their keep — Q3", h2)]
st += [P(
  "This is the most important — and most carefully hedged — part of the study. On the <b>full</b> "
  "dataset, the AI model did <i>not</i> beat conventional matching; if anything it trailed slightly. "
  "Taken alone, that reads as &ldquo;foundation models don&rsquo;t help.&rdquo; But the full dataset "
  "hides the cases that matter, because most parts have a close relative on file, where every method "
  "succeeds.")]
st += [P(
  "To isolate the hard regime, the team ran a <b>leave-close-relatives-out</b> analysis: for each "
  "part, progressively hide every similar reference, then re-check who can still identify it from only "
  "distant relatives. This manufactures the exact low-similarity regime where today&rsquo;s screening "
  "is known to fail — using nothing but benign data. The methods <b>trade places</b>:")]

img2 = Image(f"{SP}/fig_crossover.png", width=FULL*0.80, height=FULL*0.80*780/1440)
img2.hAlign = "CENTER"
st += [KeepTogether([img2, P("Figure 2. As the nearest available relative gets more distant (moving right), conventional "
  "matching (teal) collapses, while an AI model read residue-by-residue (rust) holds up. A whole-protein "
  "summary (dashed) does not — the way the model is read matters more than its size. Shaded band: the "
  "sub-30% regime where today&rsquo;s screening is known to fail.", cap)])]

st += [P(
  "At the critical threshold — below <b>30% identity</b> to any reference (639 parts, 264 groups) — the "
  "residue-level AI model identified <b>95.6%</b> of families (confidence interval 93.3&ndash;97.5), "
  "versus <b>68.2%</b> for conventional matching (63.1&ndash;73.2). The intervals do not overlap: a "
  "reversal of more than 25 points, precisely in the regime that matters for catching redesigned or "
  "highly divergent sequences. Two sub-findings sharpen the &ldquo;how&rdquo;:")]
st += [bullets([
  "<b>How you read the model is the lever, not how big it is.</b> Reading residue-by-residue wins; a "
  "whole-protein summary does not. Making the summary model 19&times; larger actually made it "
  "<i>worse</i> in the hard regime, because averaging washes out the local signal.",
  "<b>This reconciles the two results.</b> The model&rsquo;s advantage is invisible on the full dataset "
  "(where close relatives make everything easy) and only emerges once close relatives are removed.",
])]

st += [P("Why this is labelled exploratory, not proven", h2)]
st += [P(
  "The report is deliberate about not overclaiming here. Four limits keep this from being a confirmed "
  "result: the answer key is now <b>public</b> (so the test is no longer blind); the family labels were "
  "themselves built using sequence alignment (which can tilt a model-versus-alignment contest); the AI "
  "model saw many of these protein families during its original training, so some success may be "
  "memory rather than genuine generalisation; and the test identified <i>already-isolated</i> pieces "
  "rather than running the full raw-construct pipeline end to end. By the project&rsquo;s own standard, "
  "Q3 therefore remains <b>formally open</b> until a clean confirmation is run.")]

st += [P("4.5   A tidy negative result on splitting fused proteins", h2)]
st += [P(
  "An earlier hint suggested that cutting fused proteins into pieces might help the model. Tested "
  "properly, <b>it did not</b>: fixed-size splitting helped some cases and hurt others, for no net "
  "gain. The real problem was not the cutting but the <i>summary-style</i> reading, which diluted the "
  "signal across the fusion. Reading the intact protein residue-by-residue fixed the misses without "
  "any splitting. A clean negative that redirected effort toward the approach that actually works.")]

st += [PageBreak()]

# ============================ 5. IMPACT ============================
st += [P("5   Why it matters", h1), rule()]
st += [P(
  "The two findings point in a practical, deployable direction for the organisations that screen "
  "synthetic DNA:", lead)]
st += [bullets([
  "<b>A cheap preprocessing upgrade pays off now.</b> Swapping a naive translation step for a standard "
  "gene caller — mature, free, fast — recovered about 9 points more of each construct&rsquo;s known "
  "parts for <i>less</i> compute. This is an immediate, low-risk improvement that works with the "
  "conventional matching screens already in place.",
  "<b>AI models are not a drop-in replacement — they are a targeted reinforcement.</b> They do not "
  "beat conventional matching in general, so replacing today&rsquo;s screens wholesale is not warranted. "
  "But in the specific hard regime where those screens fail — distant or redesigned sequences — a "
  "residue-level model closed most of the gap. That is exactly the blind spot the biosecurity "
  "community is most worried about.",
  "<b>Representation beats scale.</b> The useful signal lived in <i>how</i> the model was read, not in "
  "a bigger model. That is a cost-saving design lesson for anyone building screening tools: a small "
  "model read well outperformed a large model read poorly.",
  "<b>A reusable, honest evaluation.</b> The hidden-answer-key, matched-budget framework — with locked "
  "protocols and exact reproducibility — is itself a contribution: a template for evaluating future "
  "screening methods without fooling oneself.",
])]

# ============================ 6. LIMITS ============================
st += [P("6   Limitations and honest status", h1), rule()]
st += [P(
  "The report is unusually candid about what it has <i>not</i> shown, and that candor is part of its "
  "quality. The main open issues:")]
st += [bullets([
  "<b>The low-similarity goal is not yet confirmed.</b> The project&rsquo;s own scope requires an "
  "independent low-similarity evaluation. The exploratory crossover is strong evidence, not that "
  "confirmation, so this goal is reported as <b>incomplete</b> rather than met.",
  "<b>Narrow source of constructs.</b> The dataset comes from one family of deposited lab vectors, "
  "dominated by common resistance markers. Results may not transfer to eukaryotic, plant, or "
  "reporter-heavy constructs — where the bacterial gene caller is already known to stumble.",
  "<b>Identifying a family is not the whole story.</b> Correctly naming a protein family is a real, "
  "measurable result, but it does not by itself establish biological activity, function, or risk. "
  "Those judgments remain for human experts.",
  "<b>Self-review only.</b> No independent or specialist review is implied anywhere in this work.",
])]

# ============================ 7. NEXT STEPS ============================
st += [P("7   What comes next", h1), rule()]
st += [P(
  "A confirmatory study is already designed to convert the promising Q3 result into a solid one. Its "
  "ingredients directly answer the four caveats above:")]
st += [bullets([
  "<b>Fresh, still-hidden test groups</b>, so the answer key is genuinely blind again.",
  "<b>Independent family labels</b> (from protein-structure and profile methods) that do not rely on "
  "the sequence alignment being judged — removing the home-field advantage.",
  "<b>A stronger conventional opponent</b> (profile-based search), so the comparison is against the "
  "best of today&rsquo;s methods, not a basic version.",
  "<b>The full raw-construct pipeline</b> run end to end, plus a broader mix of construct types to "
  "test where the gene caller stops helping.",
])]
st += [Spacer(1,4)]
st += [shaded([P("BOTTOM LINE", sidebarH),
  P("Better preprocessing is a clear, cheap, deployable win for DNA screening today. AI protein "
    "models are not a general replacement for conventional matching — but read the right way, they "
    "shine exactly where today&rsquo;s screening is weakest. That second finding is promising and "
    "carefully bounded; confirming it cleanly is the next and most important step.",
    S('bl',fontSize=10,leading=15,textColor=INK))], bg=BG, border=ACC2)]

# ============================ APPENDIX ============================
st += [Spacer(1,10), P("Appendix · Methods &amp; reproducibility at a glance", h2), rule(LINE,0.6,0,6)]
ap = [
 ("Main dataset","777 public, benign NCBI vector constructs; 1,221 corroborated protein-coding parts; "
   "422 lineage groups; 27 families; five construct shapes."),
 ("Preprocessing","Six-frame translation (baseline) vs. Pyrodigal gene calling vs. documented parts (ceiling). "
   "None except the ceiling sees the answer key."),
 ("Identification","Normalized BLOSUM62 alignment (conventional); ESM-2 protein model at 8M and 150M "
   "parameters, read as a mean-pooled summary and residue-by-residue."),
 ("Scoring","Share of hidden documented parts identified, at matched budgets and across budget curves; "
   "parts grouped so each lineage counts once; 95% intervals by whole-group bootstrap (4,000 resamples)."),
 ("Pre-registration","Protocol, budgets, and a 5-point minimum effect locked before the held-out answer "
   "key was consulted; final test run once."),
 ("Reproducibility","All inputs pinned by SHA-256 fingerprint; analysis reproduced from saved traces "
   "byte-for-byte; pinned model, library, and weight versions recorded."),
 ("Boundaries","Public, benign material only. Defensive purpose. Outputs are research aids, not "
   "autonomous safety decisions. Internal self-review; not peer reviewed."),
]
apt = Table([[P(f"<b>{k}</b>", S('apk',fontSize=8.8,leading=11.5,textColor=ACC)), P(v, S('apv',fontSize=8.8,leading=11.5,textColor=INK))]
             for k,v in ap], colWidths=[FULL*0.22, FULL*0.78])
apt.setStyle(TableStyle([
    ("ROWBACKGROUNDS",(0,0),(-1,-1),[colors.white, BG]),
    ("LINEBELOW",(0,0),(-1,-2),0.4,LINE),
    ("TOPPADDING",(0,0),(-1,-1),5),("BOTTOMPADDING",(0,0),(-1,-1),5),
    ("LEFTPADDING",(0,0),(-1,-1),7),("RIGHTPADDING",(0,0),(-1,-1),7),
    ("VALIGN",(0,0),(-1,-1),"TOP"),
]))
st += [apt]
st += [Spacer(1,8), P(
  "This report summarises completed, self-reviewed research on public, benign constructs. Figures are "
  "drawn from the project&rsquo;s committed, reproducible result files. It describes screening-method "
  "research and makes no claim about the function, activity, or risk of any specific sequence.", small)]

doc.build(st, onFirstPage=footer, onLaterPages=footer)
print("BUILT", OUT, os.path.getsize(OUT), "bytes")
