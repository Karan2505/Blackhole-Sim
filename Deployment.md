# 🚀 Deployment Guide for Karan2505

This guide will help you deploy the Gargantua Black Hole Simulator to GitHub Pages.

---

## 📋 Prerequisites

- Git installed on your computer
- GitHub account (Karan2505)
- Terminal/Command Prompt access

---

## 🛠️ Step-by-Step Deployment

### Step 1: Create the GitHub Repository

1. Go to **[github.com/new](https://github.com/new)**
2. Fill in the details:
   - **Repository name**: `gargantua`
   - **Description**: `🌌 Interstellar-quality black hole simulation with real physics`
   - **Visibility**: Select **Public** (required for free GitHub Pages)
   - ⚠️ **DO NOT** check "Add a README file" (we have our own)
   - ⚠️ **DO NOT** check "Add .gitignore" (we have our own)
3. Click **"Create repository"**

---

### Step 2: Initialize and Push Your Code

Open your terminal/command prompt in the project folder and run these commands **one by one**:

```bash
# 1. Initialize git repository
git init

# 2. Add all files to staging
git add .

# 3. Create the first commit
git commit -m "🌌 Initial release: Gargantua Black Hole Simulator v1.0.0"

# 4. Rename branch to main (if needed)
git branch -M main

# 5. Add your GitHub repository as remote
git remote add origin https://github.com/Karan2505/gargantua.git

# 6. Push to GitHub
git push -u origin main
```

If prompted for credentials:
- **Username**: Karan2505
- **Password**: Use a [Personal Access Token](https://github.com/settings/tokens) (not your password)

---

### Step 3: Enable GitHub Pages

1. Go to your repository: **[github.com/Karan2505/gargantua](https://github.com/Karan2505/gargantua)**
2. Click the **"Settings"** tab (top menu)
3. In the left sidebar, click **"Pages"**
4. Under **"Build and deployment"**:
   - **Source**: Select **"GitHub Actions"**
5. That's it! The workflow will run automatically.

---

### Step 4: Wait for Deployment

1. Go to the **"Actions"** tab in your repository
2. You'll see a workflow running called "Deploy to GitHub Pages"
3. Wait for it to complete (green checkmark ✅)
4. This usually takes 1-2 minutes

---

### Step 5: Access Your Live Site! 🎉

Your simulation is now live at:

## **🌌 [https://karan2505.github.io/gargantua/](https://karan2505.github.io/gargantua/)**

---

## 🔧 Troubleshooting

### Problem: "Pages not found" or 404 error
**Solution**: Wait 5-10 minutes. GitHub Pages can take time to propagate.

### Problem: Workflow failed
**Solution**: 
1. Go to Actions tab
2. Click on the failed workflow
3. Check the error message
4. Usually it's a permission issue - ensure Pages is set to "GitHub Actions"

### Problem: Can't push - authentication failed
**Solution**: Create a Personal Access Token:
1. Go to [github.com/settings/tokens](https://github.com/settings/tokens)
2. Click "Generate new token (classic)"
3. Select scopes: `repo` (full control)
4. Copy the token and use it as your password

### Problem: Remote origin already exists
**Solution**: 
```bash
git remote remove origin
git remote add origin https://github.com/Karan2505/gargantua.git
git push -u origin main
```

---

## 📁 What Gets Deployed

| File/Folder | Purpose |
|-------------|---------|
| `index.html` | Main interactive visualization |
| `README.md` | Project documentation |
| `LICENSE` | MIT License |
| `core/` | C++ physics engine code |
| `cuda/` | GPU acceleration kernels |
| `python/` | Python bindings |
| `docs/` | Additional documentation |
| `examples/` | Example notebooks and scripts |

---

## 🔄 Updating Your Site

To make changes and update the live site:

```bash
# Make your changes to the code, then:
git add .
git commit -m "✨ Description of your changes"
git push
```

The site will automatically redeploy in 1-2 minutes!

---

## 📱 Sharing Your Project

Share these links:

- **Live Demo**: `https://karan2505.github.io/gargantua/`
- **Source Code**: `https://github.com/Karan2505/gargantua`
- **Clone Command**: `git clone https://github.com/Karan2505/gargantua.git`

---

## 🏆 Badge for Your Profile

Add this to other README files or your GitHub profile:

```markdown
[![Gargantua](https://img.shields.io/badge/🌌_Gargantua-Black_Hole_Simulator-orange?style=for-the-badge)](https://karan2505.github.io/gargantua/)
```

Result: [![Gargantua](https://img.shields.io/badge/🌌_Gargantua-Black_Hole_Simulator-orange?style=for-the-badge)](https://karan2505.github.io/gargantua/)

---

## ✅ Deployment Checklist

- [ ] Created repository on GitHub
- [ ] Ran `git init`
- [ ] Ran `git add .`
- [ ] Ran `git commit -m "..."`
- [ ] Ran `git remote add origin ...`
- [ ] Ran `git push -u origin main`
- [ ] Enabled GitHub Pages (Settings → Pages → GitHub Actions)
- [ ] Waited for Actions workflow to complete
- [ ] Tested live site: https://karan2505.github.io/gargantua/

---

**Congratulations! 🎉 Your Interstellar Black Hole Simulator is now live!**