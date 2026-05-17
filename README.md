# 初心者向け Agent AI チュートリアル

Agent AI (Claude Code, OpenAI Codex, Gemini CLI)は何ができる？

- コンピュータ内ファイルの読み取りと編集
  - ファイルをチャット欄に添付する必要がない
  - AIの編集案を手動で反映する必要がない
- コマンド実行
  - ターミナルでできることは何でもAIにさせられる．PC内データの全消去も．

まずはインストール：
```bash
brew install gemini-cli # Gemini CLIは無料で使える
brew install claude codex # Claude Code, Codexは有料プラン必要．
```
このチュートリアルは特定のAgentに限らない内容です．

起動：
```bash
cd /path/to/agentai-tutorial
gemini
```

聞いてみよう：
> このプロジェクトのスライドファイルをブラウザで表示する方法は？

AIがファイルを読み，方法を教えてくれるはず．
