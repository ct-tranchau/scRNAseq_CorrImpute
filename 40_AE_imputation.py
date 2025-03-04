import argparse
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data
import pandas as pd
import numpy as np

# Argument parser for command-line execution
parser = argparse.ArgumentParser(description="Autoencoder for single-cell RNA-seq imputation")
parser.add_argument("--input_path", type=str, required=True, help="Path to input CSV file")
parser.add_argument("--output_path", type=str, required=True, help="Path to save imputed output CSV file")
parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs")
parser.add_argument("--batch_size", type=int, default=64, help="Batch size for training")
parser.add_argument("--lr", type=float, default=0.0001, help="Learning rate")
args = parser.parse_args()

# Step 1: Load and preprocess the single-cell RNA-seq data
print(f"Loading data from {args.input_path}...")
original_data = pd.read_csv(args.input_path, index_col=0)
transformed_data = np.log1p(original_data)  # Apply log1p transformation

# Convert data to PyTorch tensor
input_data = torch.tensor(transformed_data.values, dtype=torch.float32)

# Step 2: Create a PyTorch DataLoader
train_dataset = data.TensorDataset(input_data)
train_loader = data.DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)

# Step 3: Define the Autoencoder model
class Autoencoder(nn.Module):
    def __init__(self, input_dim, hidden_dim=400, latent_dim=200):
        super(Autoencoder, self).__init__()

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, latent_dim),
            nn.LeakyReLU(0.2)
        )
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, input_dim),
            nn.ReLU()
        )
     
    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

# Step 4: Define the loss function
def loss_function(x, x_hat):
    loss = nn.functional.mse_loss(x_hat, x, reduction='mean')
    return loss

# Step 5: Initialize the model, optimizer, and device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
input_dim = input_data.shape[1]
model = Autoencoder(input_dim=input_dim).to(device)
optimizer = optim.Adam(model.parameters(), lr=args.lr)

# Step 6: Train the model
def train(model, optimizer, train_loader, epochs, device):
    model.train()
    for epoch in range(epochs):
        overall_loss = 0
        for batch_idx, (x,) in enumerate(train_loader):
            x = x.to(device)
            optimizer.zero_grad()
            x_hat = model(x)
            loss = loss_function(x, x_hat)
            loss.backward()
            optimizer.step()
            overall_loss += loss.item()
        
        # Compute and print average loss per epoch
        average_loss = overall_loss / len(train_loader.dataset)
        print(f"Epoch {epoch + 1}/{epochs}, Average Loss: {average_loss:.4f}")

# Train the model
print("Training Autoencoder...")
train(model, optimizer, train_loader, args.epochs, device)

# Step 7: Reconstruct the data
print("Reconstructing data...")
model.eval()
with torch.no_grad():
    reconstructed_data = model(input_data.to(device))

# Inverse transformation
reconstructed_data = np.expm1(reconstructed_data.cpu().numpy())
reconstructed_data = np.round(reconstructed_data).astype(int)
reconstructed_df = pd.DataFrame(reconstructed_data, columns=original_data.columns, index=original_data.index)

# Save output
print(f"Saving imputed data to {args.output_path}...")
reconstructed_df.to_csv(args.output_path)
print("Done!")
