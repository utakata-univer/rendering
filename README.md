# ray_trace

### 動作確認

- Windows11 WSL2
- プラットフォーム x64
- g++ (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0

### コンパイル・実行

```
export CPATH=vectormath/include/
g++ -o out rayt.cpp
./out
```

### 参考

https://zenn.dev/mebiusbox/books/8d9c42883df9f6/viewer/b85221

### 備考

- 拡散反射のみを考慮したシンプルなレイトレースを原理理解のために実装.
- 影を表示したいため三角形のz負方向に大きな球面を設置.
