import argparse
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data as data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Argument parser for command-line execution
parser = argparse.ArgumentParser(description="GAN for single-cell RNA-seq imputation")
parser.add_argument("--input_path", type=str, required=True, help="Path to input CSV file")
parser.add_argument("--output_path", type=str, required=True, help="Path to save imputed output CSV file")
parser.add_argument("--loss_plot_path", type=str, required=True, help="Path to save loss plot PDF file")
parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs")
parser.add_argument("--batch_size", type=int, default=64, help="Batch size for training")
parser.add_argument("--lr_g", type=float, default=0.0001, help="Learning rate for Generator")
parser.add_argument("--lr_d", type=float, default=0.00001, help="Learning rate for Discriminator")
args = parser.parse_args()

# Step 1: Load and preprocess the single-cell RNA-seq data
print(f"Loading data from {args.input_path}...")
original_data = pd.read_csv(args.input_path, index_col=0)

# Apply log1p transformation
transformed_data = np.log1p(original_data)

# Convert the data to a PyTorch tensor
input_data = torch.tensor(transformed_data.values, dtype=torch.float32)

# Step 2: Create a PyTorch DataLoader
train_dataset = data.TensorDataset(input_data)
train_loader = data.DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)

# Step 3: Define the GAN architecture
class Generator(nn.Module):
    def __init__(self, input_dim, hidden_dim=400, latent_dim=200):
        super(Generator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, latent_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(latent_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, input_dim),
            nn.ReLU()  # non-negative reconstructed expression levels
        )

    def forward(self, x):
        return self.model(x)

class Discriminator(nn.Module):
    def __init__(self, input_dim, hidden_dim=400):
        super(Discriminator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_dim // 2, 1),
            nn.Sigmoid()  # Binary classification (real vs. fake)
        )

    def forward(self, x):
        return self.model(x)

# Step 4: Initialize the models, optimizer, and device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
input_dim = input_data.shape[1]

generator = Generator(input_dim).to(device)
discriminator = Discriminator(input_dim).to(device)

optimizer_G = optim.Adam(generator.parameters(), lr=args.lr_g)
optimizer_D = optim.Adam(discriminator.parameters(), lr=args.lr_d)

criterion_G = nn.MSELoss()  # Generator wants reconstructed data to be close to real data
criterion_D = nn.BCELoss()  # Discriminator classifies real vs. fake samples

# Step 5: Train the GAN with loss tracking
def train_gan(generator, discriminator, train_loader, epochs, device):
    generator.train()
    discriminator.train()
    
    g_losses = []
    d_losses = []

    for epoch in range(epochs):
        g_loss_total = 0
        d_loss_total = 0

        for batch_idx, (x,) in enumerate(train_loader):
            x = x.to(device)

            # Real labels = 1, Fake labels = 0
            real_labels = torch.ones((x.size(0), 1), device=device)
            fake_labels = torch.zeros((x.size(0), 1), device=device)

            # Step 1: Train the Generator
            optimizer_G.zero_grad()
            generated_data = generator(x)
            g_loss = criterion_G(generated_data, x)  # Minimize reconstruction error
            g_loss.backward()
            optimizer_G.step()
            g_loss_total += g_loss.item()

            # Step 2: Train the Discriminator
            optimizer_D.zero_grad()
            real_pred = discriminator(x)
            fake_pred = discriminator(generated_data.detach())

            d_loss_real = criterion_D(real_pred, real_labels)
            d_loss_fake = criterion_D(fake_pred, fake_labels)
            d_loss = (d_loss_real + d_loss_fake) / 2
            d_loss.backward()
            optimizer_D.step()
            d_loss_total += d_loss.item()

        avg_g_loss = g_loss_total / len(train_loader.dataset)
        avg_d_loss = d_loss_total / len(train_loader.dataset)
        g_losses.append(avg_g_loss)
        d_losses.append(avg_d_loss)

        print(f"Epoch {epoch+1}/{epochs}, Generator Loss: {avg_g_loss:.4f}, Discriminator Loss: {avg_d_loss:.4f}")

    return generator, g_losses, d_losses

# Train the GAN
print("Training GAN...")
generator, g_losses, d_losses = train_gan(generator, discriminator, train_loader, args.epochs, device)

# Step 6: Plot and save Generator and Discriminator Loss
plt.figure(figsize=(8, 5))
plt.plot(range(1, args.epochs + 1), g_losses, label="Generator Loss", marker='o')
plt.plot(range(1, args.epochs + 1), d_losses, label="Discriminator Loss", marker='s')
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.title("Generator and Discriminator Loss During Training")
plt.legend()
plt.grid()
plt.savefig(args.loss_plot_path, format="pdf")
plt.close()
print(f"Loss plot saved to {args.loss_plot_path}")

# Step 7: Generate imputed data
print("Generating imputed data...")
generator.eval()
with torch.no_grad():
    imputed_data = generator(input_data.to(device))

# Convert back to original scale
imputed_data = np.expm1(imputed_data.cpu().numpy())
imputed_data = np.round(imputed_data).astype(int)

imputed_df = pd.DataFrame(imputed_data, columns=original_data.columns, index=original_data.index)

# Step 8: Save the imputed DataFrame
print(f"Saving imputed data to {args.output_path}...")
imputed_df.to_csv(args.output_path)
print("Done!")
