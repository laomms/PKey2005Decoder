#pragma once
// =====================================================================
// v2 P-core affinity  (ported from the x86 MiniValidator that hit 0.3s)
//
// Rationale: on a hybrid CPU (e.g. i7-13700KF = 8 P-cores + 8 E-cores)
// the MITM's equal-size chunking is badly balanced -- fast P-cores finish
// and idle while slow E-core chunks lag, and total time is set by the
// slowest chunk. Restricting the worker pool to the P-cores makes every
// worker identical, so equal chunks become correctly balanced, removes
// hyperthread contention (one thread per physical P-core), and lets the
// scheduler keep them on the fast cores. E-cores are deliberately unused;
// on this latency-bound workload their stragglers cost more than their
// throughput adds.
// =====================================================================

#ifdef _WIN32
#include <windows.h>
#include <vector>

// One affinity mask per PHYSICAL performance core (highest EfficiencyClass).
// Empty if detection fails or the CPU is not hybrid (caller falls back).
inline std::vector<DWORD_PTR> GetPerformanceCoreMasks() {
    std::vector<DWORD_PTR> pCoreMasks;
    DWORD len = 0;
    if (!GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &len)
        && GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
        std::vector<BYTE> buf(len);
        if (GetLogicalProcessorInformationEx(RelationProcessorCore,
                (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)buf.data(), &len)) {
            DWORD offset = 0;
            UCHAR maxEfficiency = 0;
            while (offset < len) {
                auto info = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buf.data() + offset);
                if (info->Relationship == RelationProcessorCore
                    && info->Processor.EfficiencyClass > maxEfficiency) {
                    maxEfficiency = info->Processor.EfficiencyClass;
                }
                offset += info->Size;
            }
            offset = 0;
            while (offset < len) {
                auto info = reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buf.data() + offset);
                if (info->Relationship == RelationProcessorCore
                    && (info->Processor.EfficiencyClass == maxEfficiency || maxEfficiency == 0)) {
                    for (WORD i = 0; i < info->Processor.GroupCount; ++i)
                        pCoreMasks.push_back(info->Processor.GroupMask[i].Mask);
                }
                offset += info->Size;
            }
        }
    }
    return pCoreMasks;
}

// Pin the CALLING thread to one P-core (by index, wrapped) and boost priority.
// No-op-safe: if masks is empty, leaves the thread unpinned.
inline void PinToPerformanceCore(const std::vector<DWORD_PTR>& masks, size_t idx) {
    if (!masks.empty()) {
        SetThreadAffinityMask(GetCurrentThread(), masks[idx % masks.size()]);
    }
    SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
}

// Count of physical P-cores (0 if unknown / not hybrid).
inline int PerformanceCoreCount() {
    return (int)GetPerformanceCoreMasks().size();
}

#else // non-Windows: affinity is a no-op, counts fall back to std::thread
#include <vector>
#include <cstdint>
typedef uint64_t DWORD_PTR; // portable stand-in so shared code type-checks
inline std::vector<DWORD_PTR> GetPerformanceCoreMasks() { return {}; }
inline void PinToPerformanceCore(const std::vector<DWORD_PTR>&, size_t) {}
inline int PerformanceCoreCount() { return 0; }
#endif
