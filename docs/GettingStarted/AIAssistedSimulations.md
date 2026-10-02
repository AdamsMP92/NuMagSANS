# AI-assisted simulations

NuMagSANS ships a machine-readable usage guide so that general-purpose AI assistants can help prepare simulations without guessing the package's specialized API, file layout, units, or rotation conventions.

The assistant supports the workflow; it does not replace the physical model. NuMagSANS remains a deterministic scattering simulator, and the user remains responsible for the physical assumptions, real-space state, units, numerical resolution, and interpretation of the result.

## Load the installed context

After installing NuMagSANS, print the guide together with facts from the current installation:

```bash
python -m NuMagSANS > numagsans-context.txt
```

Attach `numagsans-context.txt` to the conversation or paste it before describing the simulation. The same information is available from Python:

```python
from NuMagSANS import llm_context

print(llm_context())
```

The installed header reports the package version, expected backend executable, supported Fourier approaches, Euler conventions, and exact output keys. This is preferable to relying on API information remembered by a general model.

For documentation-aware tools, the published context is available as:

- [`llms.txt`](https://adamsmp92.github.io/NuMagSANS/llms.txt): short index of the authoritative documentation and runnable examples.
- [`llms-full.txt`](https://adamsmp92.github.io/NuMagSANS/llms-full.txt): self-contained NuMagSANS usage guide for an AI assistant.

The locally generated context is preferred because it reflects the installed package.

## Recommended workflow

Use an assistant in five explicit stages:

1. Describe the physical system and intended observable.
2. Ask the assistant to state its assumptions and identify missing scientific choices.
3. Generate real-space data and a Python script using `NuMagSANS.write_config(...)`.
4. Review the files, units, coordinate system, rotations, q grid, polarization, scattering volume, and output selection.
5. Run a reduced test case before scaling to the full simulation.

The boundary between user and assistant should remain clear:

| The assistant can help with | The user must decide or verify |
| --- | --- |
| Directory and filename conventions | Whether the real-space state is physically appropriate |
| Configuration syntax and exact output keys | Units and physical scaling |
| StructData and RotData generation | Scattering volume and polarization geometry |
| Loop organization and result paths | q range, angular resolution, and convergence |
| Preflight consistency checks | Interpretation and scientific conclusions |

## Prompt template

The following prompt makes the expected workflow explicit:

```text
Use the attached NuMagSANS context as the authoritative API guide.

I want to simulate: <physical system and purpose>.
My real-space data represent: <atomistic moments or micromagnetic cells>.
The coordinate unit is: <unit>.
The neutron beam and polarization setup is: <description>.
The observables I need are: <channels>.

Before writing files or code:
1. summarize the proposed NuMagSANS data layers and loop structure;
2. list every physical or numerical assumption;
3. ask about any scientifically important missing value;
4. do not invent configuration keys or output names.

Then produce a minimal runnable Python script using NuMagSANS.write_config().
Include a reduced preflight run before the full calculation.
```

## Example request: rotation around the beam axis

NuMagSANS defines the neutron beam direction as the x-axis and the detector as the yz-plane. A concise request for an angular sweep is:

```text
I have one magnetic object in RealSpaceData/MagData/Object_1/m_1.csv.
Generate 36 object orientations around the neutron-beam x-axis without
duplicating the MagData file. Use an xyz RotData convention, a RotData loop,
and radians in the rotation files. Request Unpolarized_2D, SpinFlip_2D, and
SpinFlip_1D output. Show the assumptions that still require my confirmation.
```

For `RotDataConvention="xyz"`, the row `phi 0 0` represents the required active x-axis rotation. The full machine-readable guide contains a complete script for this pattern.

## Reliability rules

An assistant working with NuMagSANS should follow these rules:

- Prefer the installed context, Python introspection, current documentation, and runnable examples over remembered knowledge.
- Use `NuMagSANS.write_config(...)`; do not invent or silently omit configuration keys.
- Query `NuMagSANS.generate_all_outputs()` rather than guessing output names.
- Keep atomistic and micromagnetic scaling assumptions explicit.
- Treat RotData angles as radians and global sample-rotation angles as degrees.
- Remember that the beam is along x, not the frequently assumed z-axis.
- Validate object counts against every active StructData and RotData table.
- Do not equate successful execution with physical validation.
- Preserve generated configurations and record the package version for reproducibility.

Analytical limits, symmetry checks, grid-convergence studies, and published reference datasets remain the appropriate validation methods for scientific results.
