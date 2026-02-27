# Exercise: Setting Up Your Playground with UV

## Goal

Initialize the **`playground/`** folder as a UV-managed Python project. This is the environment you will use for **all exercises in this course**. Exercises are organized inside it by lecture and exercise number.

## Steps

### 1. Navigate to the playground folder

The `playground/` folder already exists at the root of this repository. Navigate into it:

> [!TIP] Terminal
> ```bash
> $ cd playground
> ```

### 2. Initialize the project with UV

Since the folder already exists, initialize inside it:

> [!TIP] Terminal
> ```bash
> $ uv init
> ```

### 3. Clean up files you don't need

`uv init` generates a few starter files that we won't use. Delete them to keep your playground tidy:

> [!TIP] Terminal
> ```bash
> $ rm hello.py README.md
> ```

> [!NOTE]
> We're keeping `.python-version` — UV uses this file to track which Python version the project needs. We'll update it in a later step with `uv python pin 3.11`.

Don't worry — the important file (`pyproject.toml`) stays. That's the one UV uses to manage your project.

> [!NOTE]
> Your `pyproject.toml` still has a line `readme = "README.md"` pointing to the file you just deleted. That is perfectly fine — it will not cause any errors. You can ignore it.

### 4. Create a virtual environment with Python 3.11

This creates an isolated folder (`.venv`) where all the Python packages you install will live. Every package you add with `uv add` from now on gets installed here — not on your system globally. Think of it as your project's private toolbox that you will keep using for all exercises going forward.

> [!TIP] Terminal
> ```bash
> $ uv venv --python 3.11
> ```

### 5. Activate the virtual environment

**macOS/Linux**:

> [!TIP] Terminal
> ```bash
> $ source .venv/bin/activate
> ```

**Windows**:

> [!TIP] Terminal
> ```powershell
> PS> .venv\Scripts\activate
> ```

### 6. Inspect the `pyproject.toml` (before adding dependencies)

Open the `pyproject.toml` file in your editor. Notice that the `dependencies` array is empty:

```toml
dependencies = []
```

Keep this in mind — you will compare it with the file after adding a package.

### 7. Create the exercise folder structure and copy the verification script

Exercises are organized by lecture, then by exercise number. Create the folder for this exercise and copy the verification script into it:

> [!NOTE]
> The backslash (`\`) in `Lecture\ 2` escapes the space in the folder name. Alternatively, you can wrap the path in quotes:
> ```bash
> cp "../Lecture 2/exercises/0_exercise/check_uv_setup.py" lecture_2/0_uv_setup/
> ```

> [!TIP] Terminal
> ```bash
> $ mkdir -p lecture_2/0_uv_setup
> $ cp ../Lecture\ 2/exercises/0_exercise/check_uv_setup.py lecture_2/0_uv_setup/
> ```

### 8. Add the required dependencies

> [!TIP] Terminal
> ```bash
> $ uv add cowsay
> ```

### 9. Inspect the `pyproject.toml` again (after adding dependencies)

Open `pyproject.toml` once more and compare it with what you saw in step 6.

You can view it in the terminal with:

> [!TIP] Terminal
> ```bash
> cat pyproject.toml
> ```
>
> Or open it in your text editor (VS Code, nano, etc.).

Answer the following questions:

1. What changed in the `dependencies` array?
2. What Python version constraint is listed under `[project]`?
3. What is the `name` of your project and where does it come from?

> **Hint**: `pyproject.toml` is the single source of truth for your project's metadata and dependencies. Every time you run `uv add` or `uv remove`, this file is updated automatically — along with the lockfile.

### 10. Verify your setup

Run the verification script from the `playground/` root so it picks up the shared virtual environment:

> [!TIP] Terminal
> ```bash
> $ uv run lecture_2/0_uv_setup/check_uv_setup.py
> ```

If everything is set up correctly, you should see all checks passing (including the installed version of `cowsay`). If any check fails, revisit the corresponding step above.

### 11. Understand `uv lock` and `uv sync`

You may have noticed that `uv add` already installed the package and updated the lockfile. So when are `uv lock` and `uv sync` actually useful? When you edit `pyproject.toml` by hand.

Try it out:

1. Open `pyproject.toml` and manually change the cowsay dependency to pin an older version.

   Your `pyproject.toml` currently has a line like:

   ```toml
   dependencies = [
       "cowsay>=6.1",
   ]
   ```

   > [!NOTE]
   > Your version number may differ — the important part is the `>=` operator, which means "this version or newer."

   Change the `>=` to `==` to pin the exact version, and set it to `5.0`:

   ```toml
   dependencies = [
       "cowsay==5.0",
   ]
   ```

2. Run `uv lock` to re-resolve dependencies according to your manual change, then `uv sync` to install them:

   > [!TIP] Terminal
   > ```bash
   > $ uv lock
   > $ uv sync
   > ```

3. Run the verification script again and check which version of `cowsay` is now installed:

   > [!TIP] Terminal
   > ```bash
   > $ uv run lecture_2/0_uv_setup/check_uv_setup.py
   > ```

4. Compare the version shown with the one from step 10. It should now reflect the older version you pinned.

> **Takeaway**: `uv add` handles locking and syncing for you automatically. Use `uv lock` and `uv sync` when you modify `pyproject.toml` by hand or when you clone a project and need to install from an existing lockfile.

## What your playground should look like now

After completing this exercise, your `playground/` folder should contain:

```
playground/
├── .venv/
├── pyproject.toml
├── uv.lock
└── lecture_2/
    └── 0_uv_setup/
        └── check_uv_setup.py
```

You will keep building on this environment in the next exercises.
