# Setting up UV

## Table of Contents

1. [What Is a Virtual Environment?](#1-what-is-a-virtual-environment)
2. [What Is UV, and Why Does It Exist?](#2-what-is-uv-and-why-does-it-exist)
3. [Installing UV on Different Platforms](#3-installing-uv-on-different-platforms)
4. [Initializing a Project with `uv init`](#4-initializing-a-project-with-uv-init)
5. [Managing Dependencies: `uv add` & `uv remove`](#5-managing-dependencies-uv-add-and-uv-remove)
6. [Syncing & Locking Environments](#6-syncing-and-locking-environments)
7. [Setting the Python Version](#7-setting-the-python-version-to-311)
8. [Summary](#8-summary-for-students)
9. [Exercise](#9-exercise)

## 1. What Is a Virtual Environment?

When you install Python on your computer, you get one global set of packages shared by every project. This causes problems fast: Project A needs version 1.0 of a library, but Project B needs version 2.0 — and they can't both be installed globally at the same time. For example, one analysis script might need an older version of a biology package like Biopython 1.70, while another needs the latest Biopython 1.83 with new features.

```
 WITHOUT virtual environments
 ────────────────────────────────────────────────

  Your computer
  ┌──────────────────────────────────────────┐
  │  Global Python                           │
  │  ┌────────────────────────────────────┐  │
  │  │  packages/                         │  │
  │  │    requests 1.0   ← Project A      │  │
  │  │    requests 2.0   ← Project B  ✗ CONFLICT!
  │  │    numpy 1.24                      │  │
  │  │    ...                             │  │
  │  └────────────────────────────────────┘  │
  │                                          │
  │  Project A ──┐                           │
  │  Project B ──┴──► same packages folder   │
  │  Project C ──┘                           │
  └──────────────────────────────────────────┘
```

A **virtual environment** solves this by giving each project its own isolated folder of packages. Think of it as a private box for each project: whatever you install inside it stays there and doesn't affect anything else on your system.

```
 WITH virtual environments
 ────────────────────────────────────────────────

  Your computer
  ┌──────────────────────────────────────────┐
  │                                          │
  │  Project A/                              │
  │  ├── .venv/          ← isolated env      │
  │  │   └── packages/                       │
  │  │       requests 1.0                    │
  │  │       numpy 1.24                      │
  │  └── main.py                             │
  │                                          │
  │  Project B/                              │
  │  ├── .venv/          ← isolated env      │
  │  │   └── packages/                       │
  │  │       requests 2.0  ✓ no conflict!    │
  │  │       pandas 2.1                      │
  │  └── main.py                             │
  │                                          │
  └──────────────────────────────────────────┘
```

In practice a virtual environment is just a directory (typically called `.venv`) inside your project. When you **activate** it, your terminal uses the Python and packages from that directory instead of the global ones.

```
 ACTIVATING an environment
 ────────────────────────────────────────────────

  Terminal
  ┌──────────────────────────────────────────┐
  │                                          │
  │  $ source .venv/bin/activate             │
  │                                          │
  │  Before activation:                      │
  │    python ──► /usr/bin/python (global)   │
  │                                          │
  │  After activation:                       │
  │    python ──► Project A/.venv/bin/python │
  │                                          │
  └──────────────────────────────────────────┘
```

**Why should you always use one?**

* **No conflicts** — each project has exactly the packages (and versions) it needs.
* **Clean system** — your global Python stays untouched.
* **Reproducibility** — anyone can recreate the same environment from your project files.

---

## 2. What Is UV, and Why Does It Exist?

**UV** is a fast, all-in-one Python package manager. It merges the capabilities of tools like `pip`, `pip-tools`, `virtualenv`, `pyenv`, `pipx`, `poetry`, and `twine` into a unified, high-performance tool.

For additional information, see official documentation [here](https://docs.astral.sh/uv/).

### It solves key problems such as:

* **Slow installs**: UV is **8–10× faster than pip without cache**, and **80–115× faster with a warm cache**. (A **cache** is a local store of previously downloaded files. A **warm cache** means UV already has the packages saved from a previous install, so it doesn't need to re-download them.)
* **Environment complexity**: Automatically manages virtual environments, dependencies, lockfiles, and even Python version installation.
* **Reproducibility**: Generates lockfiles (`uv.lock`) for consistent environments across machines. (A **lockfile** is a snapshot of the exact package versions installed in your project. It ensures anyone who clones your code gets identical dependencies.)
* **Ease of use**: Combines many separate tools into one simple interface.

---

## 3. Installing UV on Different Platforms

You can install UV in several ways depending on your system:

### Linux/macOS

> [!TIP] Terminal
> ```bash
> curl -LsSf https://astral.sh/uv/install.sh | sh
> #pip install uv
> #pipx install uv
> #brew install uv
> ```

### Windows

> [!TIP] Terminal
> ```powershell
> powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
> #pip install uv
> #pipx install uv
> #winget install --id=astral-sh.uv -e
> #scoop install main/uv
> ```

**After installation**, confirm it's working:

> [!NOTE]
> In terminal examples, the `$` represents your command prompt — don't type it. Only type what comes after.

> [!TIP] Terminal
> ```bash
> $ uv
> An extremely fast Python package manager.
>
> Usage: uv [OPTIONS] <COMMAND>
> ...
> ```

---

## 4. Initializing a Project with `uv init`

Before initializing a project, navigate to the directory where you want to create it. For this course, you'll use the `playground/` folder:

> [!TIP] Terminal
> ```bash
> cd path/to/EIPP_2026_LLMs_carli/playground
> pwd  # verify you're in the right place
> ```

Note that `uv init` does **not** create a virtual environment. To create one immediately after init, run `uv venv --python 3.11` (shown below).

To create a new Python project scaffold:

> [!TIP] Terminal
> ```bash
> $ uv init hello-world
> $ cd hello-world
> ```

Or, in an existing directory:

> [!TIP] Terminal
> ```bash
> $ mkdir hello-world
> $ cd hello-world
> $ uv init
> ```

**This generates**:

* `pyproject.toml`
* `README.md`, `hello.py` (with a simple hello-world script)
* `.python-version` (chooses a default Python version)
* `.gitignore` (may not be created if you are already inside a git repository)

To create the virtual environment, run:

> [!TIP] Terminal
> ```bash
> $ uv venv --python 3.11
> ```

This creates a `.venv` directory in your project with Python 3.11. If that version isn't installed on your system, uv will download it automatically.

You can then activate it:

* **macOS/Linux**:

  > [!TIP] Terminal
  > ```bash
  > $ source .venv/bin/activate
  > ```

* **Windows**:

  > [!TIP] Terminal
  > ```powershell
  > PS> .venv\Scripts\activate
  > ```

---

## 5. Managing Dependencies: `uv add` & `uv remove`

* **Add a dependency**:

  > [!TIP] Terminal
  > ```bash
  > $ uv add requests
  > ```

  Optionally specify versions:

  > [!TIP] Terminal
  > ```bash
  > $ uv add "requests==2.31.0"
  > ```

  Or add from Git:

  > [!TIP] Terminal
  > ```bash
  > $ uv add git+https://github.com/psf/requests
  > ```

  For `requirements.txt`:

  > [!TIP] Terminal
  > ```bash
  > $ uv add -r requirements.txt -c constraints.txt
  > ```

* **Remove a dependency**:

  > [!TIP] Terminal
  > ```bash
  > $ uv remove requests
  > ```

These commands update your `pyproject.toml`, the lockfile, and the environment (unless overridden).

> [!NOTE]
> **Spaces in file paths:** If a folder name contains a space (like `Lecture 2`), you need to either escape the space with a backslash (`Lecture\ 2`) or wrap the path in quotes (`"Lecture 2"`). Both work the same way.

---

## 6. Syncing & Locking Environments

* **Lock dependencies**:

  > [!TIP] Terminal
  > ```bash
  > $ uv lock
  > ```

* **Sync your environment**:

  > [!TIP] Terminal
  > ```bash
  > $ uv sync
  > ```

This ensures your `.venv` matches the lockfile. If you haven't activated the virtual environment yet, see the activation instructions in [Section 4](#4-initializing-a-project-with-uv-init).

> Tip from the community:
> "If you clone an existing repo, this environment can be rebuilt with `uv sync`. uv auto‑syncs on commands like `uv run`."

---

## 7. Setting the Python Version

To pin your project to Python 3.11:

> [!TIP] Terminal
> ```bash
> $ uv python pin 3.11
> ```

This updates `.python-version` accordingly. If you already created your virtual environment with `uv venv --python 3.11` (as shown in [Section 4](#4-initializing-a-project-with-uv-init)), the pin and venv are already aligned.

---

## 8. Summary for Students

1. **What & Why**: UV is a fast, all-in-one Python project manager that unifies environments, dependencies, scripts, and version management.
2. **Install** via script, Python, pipx, or brew depending on your system.
3. **Initialize** your project using `uv init`, then immediately create a virtual environment with `uv venv --python 3.11`.
4. **Add/remove** dependencies easily (`uv add` / `uv remove`).
5. **Lock** and **sync** your environment for reproducibility.
6. **Pin Python version** to 3.11 for consistency.

## 9. Exercise

Time to put this into practice! You will set up the **`playground/`** folder at the root of this repository as your working Python project. This is the single environment you will use for **all exercises across all lectures** — organized by lecture and exercise number.

See [exercise_uv.md](exercises/0_exercise/exercise_uv.md) for the hands-on exercise.
