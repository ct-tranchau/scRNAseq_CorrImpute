import argparse
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Argument parser for command-line execution
parser = argparse.ArgumentParser(description="VAE for single-cell RNA-seq imputation")
parser.add_argument("--input_path", type=str, required=True, help="Path to input CSV file")
parser.add_argument("--output_path", type=str, required=True, help="Path to save imputed output CSV file")
parser.add_argument("--loss_plot_path", type=str, required=True, help="Path to save loss plot PDF file")
parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs")
parser.add_argument("--batch_size", type=int, default=64, help="Batch size for training")
parser.add_argument("--lr", type=float, default=0.0001, help="Learning rate")
args = parser.parse_args()

# Step 1: Load and preprocess the single-cell RNA-seq data
print(f"Loading data from {args.input_path}...")
original_data = pd.read_csv(args.input_path, index_col=0)
transformed_data = np.log1p(original_data)  # Apply log1p transformation
input_data = torch.tensor(transformed_data.values, dtype=torch.float32)

# Step 2: Create a PyTorch DataLoader
train_dataset = data.TensorDataset(input_data)
train_loader = data.DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)

# Step 3: Define the VAE model
class scVAE(nn.Module):
    def __init__(self, input_dim, hidden_dim=400, latent_dim=200, device='cpu'):
        super(scVAE, self).__init__()
        self.device = device

        # Encoder
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, latent_dim),
            nn.LeakyReLU(0.2)
        )
        
        self.mean_layer = nn.Linear(latent_dim, latent_dim)
        self.logvar_layer = nn.Linear(latent_dim, latent_dim)
        
        # Decoder
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, input_dim),
            nn.ReLU()
        )
     
    def encode(self, x):
        x = self.encoder(x)
        mean, logvar = self.mean_layer(x), self.logvar_layer(x)
        return mean, logvar

    def reparameterization(self, mean, logvar):
        std = torch.exp(0.5 * logvar)
        epsilon = torch.randn_like(std).to(self.device)
        z = mean + std * epsilon
        return z

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        mean, logvar = self.encode(x)
        z = self.reparameterization(mean, logvar)
        x_hat = self.decode(z)
        return x_hat, mean, logvar

# Step 4: Loss function
def loss_function(x, x_hat, mean, log_var):
    reproduction_loss = nn.functional.mse_loss(x_hat, x, reduction='none').sum(dim=-1).mean()
    epsilon = 1e-6
    KLD = -0.5 * torch.sum(1 + log_var - mean.pow(2) - (log_var + epsilon).exp(), dim=-1).mean()
    return reproduction_loss + KLD

# Step 5: Initialize the model, optimizer, and device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
input_dim = input_data.shape[1]
model = scVAE(input_dim=input_dim, device=device).to(device)
optimizer = optim.Adam(model.parameters(), lr=args.lr)

# Step 6: Train the model
def train(model, optimizer, train_loader, epochs, device):
    model.train()
    loss_history = []

    for epoch in range(epochs):
        overall_loss = 0
        for batch_idx, (x,) in enumerate(train_loader):
            x = x.to(device)
            optimizer.zero_grad()
            x_hat, mean, log_var = model(x)
            loss = loss_function(x, x_hat, mean, log_var)
            loss.backward()
            optimizer.step()
            overall_loss += loss.item()
        
        average_loss = overall_loss / len(train_loader.dataset)
        loss_history.append(average_loss)
        print(f"Epoch {epoch + 1}/{epochs}, Average Loss: {average_loss:.4f}")

    return loss_history

# Train the model
print("Training VAE...")
loss_history = train(model, optimizer, train_loader, args.epochs, device)

# Step 7: Plot and save loss
plt.figure(figsize=(8, 5))
plt.plot(range(1, args.epochs + 1), loss_history, marker='o', linestyle='-')
plt.xlabel("Epochs")
plt.ylabel("Loss")
plt.title("Training Loss over Epochs")
plt.grid(True)
plt.savefig(args.loss_plot_path)
print(f"Loss plot saved to {args.loss_plot_path}")

# Step 8: Reconstruct the data
print("Reconstructing data...")
model.eval()
with torch.no_grad():
    reconstructed_data, _, _ = model(input_data.to(device))

# Inverse transformation
reconstructed_data = np.expm1(reconstructed_data.cpu().numpy())
reconstructed_data = np.round(reconstructed_data).astype(int)
reconstructed_df = pd.DataFrame(reconstructed_data, columns=original_data.columns, index=original_data.index)

# Save output
print(f"Saving imputed data to {args.output_path}...")
reconstructed_df.to_csv(args.output_path)
print("Done!")
