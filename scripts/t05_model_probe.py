import torch
from pathlib import Path
from transformers import AutoTokenizer, AutoModel

model_id = "facebook/esm2_t6_8M_UR50D"

# Read a local, benign protein FASTA record
ROOT = Path(__file__).resolve().parents[1]
fasta_path = ROOT / "data" / "raw" / "P42212.fasta"

lines = fasta_path.read_text(encoding="utf-8").splitlines()
sequence = "".join(
    line.strip() for line in lines if not line.startswith(">")
)

tokenizer = AutoTokenizer.from_pretrained(model_id)
model = AutoModel.from_pretrained(model_id)
model.eval()

inputs = tokenizer(sequence, return_tensors="pt")

with torch.inference_mode():
    outputs = model(**inputs)

print("Amino acids:", len(sequence))
print("Output shape:", outputs.last_hidden_state.shape)
