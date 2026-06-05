# AGENTS.md — Coding Guidelines for Dinara

## Coding Discipline (Karpathy Rules)

### 1. Think Before Coding
Don't assume. Don't hide confusion. Surface tradeoffs.
- State assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them — don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.

### 2. Simplicity First
Minimum code that solves the problem. Nothing speculative.
- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- If you write 200 lines and it could be 50, rewrite it.

### 3. Surgical Changes
Touch only what you must. Clean up only your own mess.
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- Remove imports/variables/functions that YOUR changes made unused.
- Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution
Define success criteria. Loop until verified.
- Transform tasks into verifiable goals with concrete checks.
- For multi-step tasks, state a brief plan with verification steps.
- Strong success criteria let you loop independently.

## Lesson Learned: Additive vs Destructive Changes

When improving SV call sizing, **never replace** an existing estimate destructively. Always emit both the original and the refined estimate, letting the downstream matcher (truvari `--pick multi`) choose the best. The V36x regression demonstrated this: destructive CIGAR refinement in indirect-covdrop replaced accurate multi-repeat-unit estimates with single-repeat-unit CIGAR sizes, causing -368 TP at pctsize=0.7. The fix (V36y) removed the destructive path and kept only the additive `cigarRefineInsCalls` lambda, recovering all losses and gaining +90 TP over baseline.
