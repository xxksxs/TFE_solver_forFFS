# 材料和单位检查

材料参数和单位体系错误会直接影响波数、截止频率和场幅度。

## 检查项

- 内部频率使用 Hz。
- 坐标解释与 NGMesh 文件一致。
- 相对介电常数、磁导率和电导率来自 AEDT 材料块。
- 有损材料产生合理衰减。

## 常见风险

- 频率显示使用 GHz，但内部计算误用 GHz。
- 网格坐标单位与波数计算不一致。
- 材料参数未正确从 AEDT 映射到体单元。
- 损耗项符号或复数系数错误。

## 修改建议

涉及材料或单位修改时，应同步检查：

- `../../include/bpfem/core/Types.hpp`
- `../../src/io/AEDTParser.cpp`
- `../../src/io/NGMeshParser.cpp`
- `../../src/fem/FEMAssembler.cpp`
