/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_DOMAIN_H_
#define HIGHS_DOMAIN_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <set>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolver.h"
#include "util/HighsCDouble.h"
#include "util/HighsRbTree.h"

class HighsCutPool;
class HighsConflictPool;
class HighsObjectiveFunction;

class HighsDomain {
 public:
  struct Reason {
    HighsInt type;
    HighsInt index;

    enum {
      kBranching = -1,
      kUnknown = -2,
      kModelRowUpper = -3,
      kModelRowLower = -4,
      kCliqueTable = -5,
      kConflictingBounds = -6,
      kObjective = -7,
    };
    static Reason branching() { return Reason{kBranching, 0}; }
    static Reason unspecified() { return Reason{kUnknown, 0}; }
    static Reason cliqueTable(HighsInt col, HighsInt val) {
      return Reason{kCliqueTable, 2 * col + val};
    }
    static Reason modelRowUpper(HighsInt row) {
      return Reason{kModelRowUpper, row};
    }
    static Reason modelRowLower(HighsInt row) {
      return Reason{kModelRowLower, row};
    }
    static Reason cut(HighsInt cutpool, HighsInt cut) {
      return Reason{cutpool, cut};
    }
    static Reason conflictingBounds(HighsInt pos) {
      return Reason{kConflictingBounds, pos};
    }
    static Reason objective() { return Reason{kObjective, 0}; }
  };

  class ConflictSet {
    friend class HighsDomain;
    HighsDomain& localdom;
    HighsDomain& globaldom;

   public:
    struct LocalDomChg {
      HighsInt pos;
      mutable HighsDomainChange domchg;

      bool operator<(const LocalDomChg& other) const { return pos < other.pos; }
    };

    ConflictSet(HighsDomain& localdom);

    void conflictAnalysis(HighsConflictPool& conflictPool);
    void conflictAnalysis(const HighsInt* proofinds, const double* proofvals,
                          HighsInt prooflen, double proofrhs,
                          HighsConflictPool& conflictPool);

   private:
    std::set<LocalDomChg> reasonSideFrontier;
    std::set<LocalDomChg> reconvergenceFrontier;
    std::vector<std::set<LocalDomChg>::iterator> resolveQueue;
    std::vector<LocalDomChg> resolvedDomainChanges;

    struct ResolveCandidate {
      double delta;
      double baseBound;
      double prio;
      HighsInt boundPos;
      HighsInt valuePos;

      bool operator<(const ResolveCandidate& other) const {
        if (prio > other.prio) return true;
        if (other.prio > prio) return false;

        return boundPos < other.boundPos;
      }
    };

    std::vector<ResolveCandidate> resolveBuffer;

    void pushQueue(std::set<LocalDomChg>::iterator domchgPos);
    std::set<LocalDomChg>::iterator popQueue();
    void clearQueue();
    HighsInt queueSize() const;
    bool resolvable(HighsInt domChgPos) const;

    HighsInt resolveDepth(std::set<LocalDomChg>& frontier, HighsInt depthLevel,
                          HighsInt stopSize, HighsInt minResolve = 0,
                          bool increaseConflictScore = false);

    HighsInt computeCuts(HighsInt depthLevel, HighsConflictPool& conflictPool);

    bool explainInfeasibility();

    bool explainInfeasibilityConflict(const HighsDomainChange* conflict,
                                      HighsInt len);

    bool explainInfeasibilityLeq(const HighsInt* inds, const double* vals,
                                 HighsInt len, double rhs, double minActivity);

    bool explainInfeasibilityGeq(const HighsInt* inds, const double* vals,
                                 HighsInt len, double rhs, double maxActivity);

    bool explainBoundChange(const std::set<LocalDomChg>& currentFrontier,
                            LocalDomChg domchg);

    // bool explainBoundChange(HighsInt pos) {
    //   return explainBoundChange(LocalDomChg{pos,
    //   localdom.domchgstack_[pos]});
    // }

    bool explainBoundChangeConflict(const LocalDomChg& domchg,
                                    const HighsDomainChange* conflict,
                                    HighsInt len);

    bool explainBoundChangeLeq(const std::set<LocalDomChg>& currentFrontier,
                               const LocalDomChg& domChg, const HighsInt* inds,
                               const double* vals, HighsInt len, double rhs,
                               double minActivity);

    bool explainBoundChangeGeq(const std::set<LocalDomChg>& currentFrontier,
                               const LocalDomChg& domChg, const HighsInt* inds,
                               const double* vals, HighsInt len, double rhs,
                               double maxActivity);

    bool resolveLinearLeq(HighsCDouble M, double Mlower, const double* vals);

    bool resolveLinearGeq(HighsCDouble M, double Mupper, const double* vals);
  };

  struct CutpoolPropagation {
    HighsInt cutpoolindex;
    HighsDomain* domain;
    HighsCutPool* cutpool;
    std::vector<HighsCDouble> activitycuts_;
    std::vector<HighsInt> activitycutsinf_;
    std::vector<uint8_t> propagatecutflags_;
    std::vector<HighsInt> propagatecutinds_;
    std::vector<double> capacityThreshold_;

    CutpoolPropagation(HighsInt cutpoolindex, HighsDomain* domain,
                       HighsCutPool& cutpool);

    CutpoolPropagation(const CutpoolPropagation& other);

    ~CutpoolPropagation();

    void recomputeCapacityThreshold(HighsInt cut);

    void cutAdded(HighsInt cut, bool propagate);

    void cutDeleted(HighsInt cut, bool deletedOnlyForPropagation = false);

    void markPropagateCut(HighsInt cut);

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound);
  };

  struct ConflictPoolPropagation {
    HighsInt conflictpoolindex;
    HighsDomain* domain;
    HighsConflictPool* conflictpool_;
    std::vector<HighsInt> colLowerWatched_;
    std::vector<HighsInt> colUpperWatched_;
    std::vector<uint8_t> conflictFlag_;
    std::vector<HighsInt> propagateConflictInds_;

    struct WatchedLiteral {
      HighsDomainChange domchg = {0.0, -1, HighsBoundType::kLower};
      HighsInt prev = -1;
      HighsInt next = -1;
    };

    std::vector<WatchedLiteral> watchedLiterals_;

    ConflictPoolPropagation(HighsInt conflictpoolindex, HighsDomain* domain,
                            HighsConflictPool& cutpool);

    ConflictPoolPropagation(const ConflictPoolPropagation& other);

    ~ConflictPoolPropagation();

    void linkWatchedLiteral(HighsInt linkPos);

    void unlinkWatchedLiteral(HighsInt linkPos);

    void conflictAdded(HighsInt conflict);

    void conflictDeleted(HighsInt conflict);

    void markPropagateConflict(HighsInt conflict);

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

    void propagateConflict(HighsInt conflict);
  };

 private:
  struct ObjectivePropagation {
    HighsDomain* domain = nullptr;
    const HighsObjectiveFunction* objFunc;
    const double* cost;
    HighsCDouble objectiveLower;
    HighsInt numInfObjLower;
    double capacityThreshold;
    bool isPropagated;

    struct ObjectiveContribution {
      double contribution;
      HighsInt col;
      HighsInt partition;
      highs::RbTreeLinks<HighsInt> links;
    };

    class ObjectiveContributionTree;

    std::vector<ObjectiveContribution> objectiveLowerContributions;
    std::vector<std::pair<HighsInt, HighsInt>> contributionPartitionSets;
    std::vector<double> propagationConsBuffer;
    struct PartitionCliqueData {
      double multiplier;
      HighsInt rhs;
      bool changed;
    };

    std::vector<PartitionCliqueData> partitionCliqueData;

    ObjectivePropagation() {
      objFunc = nullptr;
      cost = nullptr;
      objectiveLower = 0.0;
      numInfObjLower = 0;
      capacityThreshold = 0.0;
      isPropagated = false;
    }
    ObjectivePropagation(HighsDomain* domain);

    bool isActive() const { return domain != nullptr; }

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

    bool shouldBePropagated() const;

    void propagate();

    void debugCheckObjectiveLower() const;

    // construct the proot constraint at the time when the domain change stack
    // had the given size
    void getPropagationConstraint(HighsInt domchgStackSize, const double*& vals,
                                  const HighsInt*& inds, HighsInt& len,
                                  double& rhs, HighsInt domchgCol = -1);

   private:
    void recomputeCapacityThreshold();
  };

  std::vector<uint8_t> changedcolsflags_;
  std::vector<HighsInt> changedcols_;

  std::vector<std::pair<HighsInt, HighsInt>> propRowNumChangedBounds_;

  std::vector<HighsDomainChange> domchgstack_;
  std::vector<Reason> domchgreason_;
  std::vector<std::pair<double, HighsInt>> prevboundval_;

  std::vector<HighsCDouble> activitymin_;
  std::vector<HighsCDouble> activitymax_;
  std::vector<HighsInt> activitymininf_;
  std::vector<HighsInt> activitymaxinf_;
  std::vector<double> capacityThreshold_;
  std::vector<uint8_t> propagateflags_;
  std::vector<HighsInt> propagateinds_;
  ObjectivePropagation objProp_;


  // row lower and upper, length = 2 * rownum
  std::vector<char> redundantPropagateflags_;
  std::vector<HighsInt> redundantPropagateinds_;
  std::vector<std::pair<HighsInt, bool>> zeroCostFixedVariables_;
  // vectors to record if a constraint is redundant in probing x_k = 1
  std::vector<HighsInt> firstOccuranceInds_;
  std::vector<char> firstOccuranceFlag_;
  // vectors to record if a constraint is redundant in probing both x_k = 1 and x_k = 0
  // If so, we may derive some GDF fixings for variables in this constraint. 
  std::vector<HighsInt> twoSideRedundantInds_;
  std::vector<char> redundantResidueflags_;
  std::vector<HighsInt> redundantResidueinds_;
  std::vector<double> redundantAbsResidue_;
  std::unordered_map<HighsInt, std::unordered_set<HighsInt>> permanentLowerLock_;
  std::unordered_map<HighsInt, std::unordered_set<HighsInt>> permanentUpperLock_;

  HighsInt probingStatusSide = 0;
  bool startZeroCostFixing;

  std::vector<HighsInt> tmpColLoLock_;
  std::vector<HighsInt> tmpColUpLock_;
  std::vector<HighsInt> involvedVars;
  std::vector<char> indsVars;
  void clearInvolved(HighsInt start) {
    for (const auto x : involvedVars)
      indsVars[x] = false;
    involvedVars.clear();
  }

  std::vector<HighsInt> varTightenInds_;
  std::vector<char> varTightenFlag_;
  std::vector<double> tightestColLower;
  std::vector<double> tightestColUpper;

  void insertVarsBchg(const HighsDomainChange& bchg) {
    HighsInt col = bchg.column;
    if (!varTightenFlag_[col]) {
      varTightenFlag_[col] = true;
      varTightenInds_.push_back(col);
    }
    if (bchg.boundtype == HighsBoundType::kLower)
      tightestColLower[col] = std::max(tightestColLower[col], bchg.boundval);
    else
      tightestColUpper[col] = std::min(tightestColUpper[col], bchg.boundval);
  }

  bool isBchgTightenType(const HighsDomainChange& bchg) {
    if (bchg.boundtype == HighsBoundType::kLower && bchg.boundval > col_lower_[bchg.column])
      return true;
    if (bchg.boundtype == HighsBoundType::kUpper && bchg.boundval < col_upper_[bchg.column])
      return true;
    return false;
  }

  HighsMipSolver* mipsolver;

 private:
  std::deque<CutpoolPropagation> cutpoolpropagation;
  std::deque<ConflictPoolPropagation> conflictPoolPropagation;

  bool infeasible_ = false;
  bool inProbing_ = false;
  bool probingPrepared_ = false;
  Reason infeasible_reason;
  HighsInt infeasible_pos;

  void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

  void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

  void updateThresholdLbChange(HighsInt col, double newbound, double val,
                               double& threshold) const;

  void updateThresholdUbChange(HighsInt col, double newbound, double val,
                               double& threshold) const;

  void recomputeCapacityThreshold(HighsInt row);

  void updateRedundantRows(HighsInt row);

  double doChangeBound(const HighsDomainChange& boundchg);

  std::vector<HighsInt> colLowerPos_;
  std::vector<HighsInt> colUpperPos_;
  std::vector<HighsInt> branchPos_;
  HighsHashTable<HighsInt> redundantRows_;
  bool recordRedundantRows_ = false;

 public:
  std::vector<double> col_lower_;
  std::vector<double> col_upper_;
  void markInProbing() { inProbing_ = 1; }
  void markNoProbing() { inProbing_ = 0; }
  HighsInt zeroCostStartPos_ = kHighsIInf32;
  bool allowDFProbing = true;
 

  HighsDomain(HighsMipSolver& mipsolver);

  HighsDomain(const HighsDomain& other)
      : changedcolsflags_(other.changedcolsflags_),
        changedcols_(other.changedcols_),
        domchgstack_(other.domchgstack_),
        domchgreason_(other.domchgreason_),
        prevboundval_(other.prevboundval_),
        activitymin_(other.activitymin_),
        activitymax_(other.activitymax_),
        activitymininf_(other.activitymininf_),
        activitymaxinf_(other.activitymaxinf_),
        capacityThreshold_(other.capacityThreshold_),
        propagateflags_(other.propagateflags_),
        propagateinds_(other.propagateinds_),
        objProp_(other.objProp_),
        mipsolver(other.mipsolver),
        cutpoolpropagation(other.cutpoolpropagation),
        conflictPoolPropagation(other.conflictPoolPropagation),
        infeasible_(other.infeasible_),
        infeasible_reason(other.infeasible_reason),
        infeasible_pos(other.infeasible_pos),
        colLowerPos_(other.colLowerPos_),
        colUpperPos_(other.colUpperPos_),
        branchPos_(other.branchPos_),
        col_lower_(other.col_lower_),
        col_upper_(other.col_upper_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.domain = this;
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.domain = this;
    if (objProp_.domain) objProp_.domain = this;
  }

  HighsDomain& operator=(const HighsDomain& other) {
    changedcolsflags_ = other.changedcolsflags_;
    changedcols_ = other.changedcols_;
    domchgstack_ = other.domchgstack_;
    domchgreason_ = other.domchgreason_;
    prevboundval_ = other.prevboundval_;
    activitymin_ = other.activitymin_;
    activitymax_ = other.activitymax_;
    activitymininf_ = other.activitymininf_;
    activitymaxinf_ = other.activitymaxinf_;
    capacityThreshold_ = other.capacityThreshold_;
    propagateflags_ = other.propagateflags_;
    propagateinds_ = other.propagateinds_;
    objProp_ = other.objProp_;
    mipsolver = other.mipsolver;
    cutpoolpropagation = other.cutpoolpropagation;
    conflictPoolPropagation = other.conflictPoolPropagation;
    infeasible_ = other.infeasible_;
    infeasible_reason = other.infeasible_reason;
    colLowerPos_ = other.colLowerPos_;
    colUpperPos_ = other.colUpperPos_;
    branchPos_ = other.branchPos_;
    col_lower_ = other.col_lower_;
    col_upper_ = other.col_upper_;
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.domain = this;
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.domain = this;
    if (objProp_.domain) objProp_.domain = this;
    return *this;
  }

  void computeMinActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const double* ARvalue, HighsInt& ninfmin,
                          HighsCDouble& activitymin) const;

  void computeMaxActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const double* ARvalue, HighsInt& ninfmax,
                          HighsCDouble& activitymax) const;

  double adjustedUb(HighsInt col, HighsCDouble boundVal, bool& accept) const;

  double adjustedLb(HighsInt col, HighsCDouble boundVal, bool& accept) const;

  HighsInt propagateRowUpper(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rupper,
                             const HighsCDouble& minactivity, HighsInt ninfmin,
                             HighsDomainChange* boundchgs) const;

  HighsInt propagateRowLower(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rlower,
                             const HighsCDouble& maxactivity, HighsInt ninfmax,
                             HighsDomainChange* boundchgs) const;

  const std::vector<HighsInt>& getChangedCols() const { return changedcols_; }

  void addCutpool(HighsCutPool& cutpool);

  void addConflictPool(HighsConflictPool& conflictPool);

  void clearChangedCols() {
    for (HighsInt i : changedcols_) changedcolsflags_[i] = 0;
    changedcols_.clear();
  }

  void removeContinuousChangedCols() {
    for (HighsInt i : changedcols_)
      changedcolsflags_[i] = mipsolver->isColIntegral(i);

    changedcols_.erase(
        std::remove_if(changedcols_.begin(), changedcols_.end(),
                       [&](HighsInt i) { return !isChangedCol(i); }),
        changedcols_.end());
  }

  void clearChangedCols(size_t start) {
    for (size_t i = start; i != changedcols_.size(); ++i)
      changedcolsflags_[changedcols_[i]] = 0;

    changedcols_.resize(start);
  }

  bool isChangedCol(HighsInt col) const { return changedcolsflags_[col] != 0; }

  void markPropagate(HighsInt row);

  // mark constraints as redundant, to compute the set of non-redundant rows
  void markRedundantPropagate(HighsInt row, bool isUpper);
  // reset redundant status due to infeasibility encountered in propagation
  void cancelRedundantPropagate(HighsInt row, bool isUpper);
  // check if the "<=" side of row-th constraint is redundant
  bool isUpperRedundant(HighsInt row);
  // check if the ">=" side of row-th constraint is redundant
  bool isLowerRedundant(HighsInt row);
  /* vectors to store fixing directions for variables with zero objective coefficients
     Possible values:
      1  (only allowed to be fixed to upper bound) 
      -1 (only allowed to be fixed to lower bound) 
      0  (not specified yet) 
  */
  int* zeroCostVarsDirection_;

  HighsInt getRedundantSize() {
    return redundantPropagateinds_.size();
  }

  // clear row redundant occurance information
  void clearOccruance() {
    for (auto x : firstOccuranceInds_)
      firstOccuranceFlag_[x] = false;
    firstOccuranceInds_.clear();
    twoSideRedundantInds_.clear();
  }

  // clear the records of the tightest bounds derived in probing
  void clearAllVariableTight() {
    for (auto x : varTightenInds_) {
      varTightenFlag_[x] = false;
      tightestColLower[x] = col_lower_[x];
      tightestColUpper[x] = col_upper_[x];
    }
    varTightenInds_.clear();
  }
  // status = -1 (no probing), 0 (probing at 0), 1(probing at 1)
  void setProbingStatusFlag(HighsInt a) {
    probingStatusSide = a;
  }

  // clear the lock (which can be removed for dual fixing) that are identified by GDF
  void clearPermanentLock() {
    permanentLowerLock_.clear();
    permanentUpperLock_.clear();
  }
  // clear redundant flags and vectors for constraints
  void clearRedundant();
  // compute redundant residues for GDF
  void frozenRedundantResidue();
  // clear redundant resiude data
  void clearRedundantResidue();
  // apply GDF on the probing variable itself
  void probingVariableGDF(HighsInt col, bool val);
  // check if we can remove locks by Theorem 1, based on probing
  void updateGDFInfo();
  // apply GDF after probing both x_k = 1 and x_k = 0
  int finalRoundDualFixing();
  // debug function
  bool cumulateReduction();

  bool isActive(const HighsDomainChange& domchg) const {
    return domchg.boundtype == HighsBoundType::kLower
               ? domchg.boundval <= col_lower_[domchg.column]
               : domchg.boundval >= col_upper_[domchg.column];
  }

  void markPropagateCut(Reason reason);

  void setupObjectivePropagation() { objProp_ = ObjectivePropagation(this); }

  void computeRowActivities();

  void initializeDualfixingUtils();

  void markInfeasible(Reason reason = Reason::unspecified()) {
    infeasible_ = true;
    infeasible_pos = domchgstack_.size();
    infeasible_reason = reason;
  }

  bool infeasible() const { return infeasible_; }

  bool inProbing() const { return inProbing_; }
  void cumulateBchg(int n);
 
  /* If a bound change tightens an infinite bound to a finite bound, we do not count it as an implication; 
      for example u_j = +\inf (global bound) and u^1_j = 0 (probing from x_k = 1),
      then the bchg reads x_j \le (+\inf) x_j, which will be rejected by HiGHS.
      See line 153 and 163 in HighsImplications.cpp
  */
  bool isFiniteBoundTightening(HighsDomainChange& boundchg);
  // increase the number of bchg by n
  void cumulateIter(int n);

  void changeBound(HighsDomainChange boundchg,
                   Reason reason = Reason::branching());

  void changeBound(HighsBoundType boundtype, HighsInt col, double boundval,
                   Reason reason = Reason::branching()) {
    changeBound({boundval, col, boundtype}, reason);
  }

  bool checkChangeBound(HighsBoundType boundtype, HighsInt col,
                        HighsCDouble boundval, Reason reason);

  void fixCol(HighsInt col, double val, Reason reason = Reason::unspecified()) {
    if (kAllowDeveloperAssert) {
      assert(infeasible_ == 0);
    }
    if (col_lower_[col] < val) {
      changeBound({val, col, HighsBoundType::kLower}, reason);
      if (infeasible_ == 0) propagate();
    }

    if (infeasible_ == 0 && col_upper_[col] > val)
      changeBound({val, col, HighsBoundType::kUpper}, reason);
  }

  void backtrackToGlobal();

  HighsDomainChange backtrack();

  const std::vector<HighsInt>& getBranchingPositions() const {
    return branchPos_;
  }

  const std::vector<std::pair<double, HighsInt>>& getPreviousBounds() const {
    return prevboundval_;
  }

  const std::vector<HighsDomainChange>& getDomainChangeStack() const {
    return domchgstack_;
  }

  const std::vector<Reason>& getDomainChangeReason() const {
    return domchgreason_;
  }

  double getObjectiveLowerBound() const {
    if (objProp_.isActive() && objProp_.numInfObjLower == 0)
      return double(objProp_.objectiveLower);

    return -kHighsInf;
  }

  void getCutoffConstraint(const double*& vals, const HighsInt*& inds,
                           HighsInt& len, double& rhs) {
    objProp_.getPropagationConstraint(domchgstack_.size(), vals, inds, len,
                                      rhs);
  }

  HighsInt getNumDomainChanges() const { return domchgstack_.size(); }

  bool colBoundsAreGlobal(HighsInt col) const {
    return colLowerPos_[col] == -1 && colUpperPos_[col] == -1;
  }

  HighsInt getBranchDepth() const { return branchPos_.size(); }

  std::vector<HighsDomainChange> getReducedDomainChangeStack(
      std::vector<HighsInt>& branchingPositions) const {
    std::vector<HighsDomainChange> reducedstack;
    reducedstack.reserve(domchgstack_.size());
    branchingPositions.reserve(branchPos_.size());
    for (HighsInt i = 0; i < (HighsInt)domchgstack_.size(); ++i) {
      // keep only the tightest bound change for each variable
      if ((domchgstack_[i].boundtype == HighsBoundType::kLower &&
           colLowerPos_[domchgstack_[i].column] != i) ||
          (domchgstack_[i].boundtype == HighsBoundType::kUpper &&
           colUpperPos_[domchgstack_[i].column] != i))
        continue;

      if (domchgreason_[i].type == Reason::kBranching)
        branchingPositions.push_back(reducedstack.size());
      else {
        HighsInt k = i;
        while (prevboundval_[k].second != -1) {
          k = prevboundval_[k].second;
          if (domchgreason_[k].type == Reason::kBranching) {
            branchingPositions.push_back(reducedstack.size());
            break;
          }
        }
      }

      reducedstack.push_back(domchgstack_[i]);
    }

    reducedstack.shrink_to_fit();
    return reducedstack;
  }

  void setDomainChangeStack(const std::vector<HighsDomainChange>& domchgstack);

  void setDomainChangeStack(const std::vector<HighsDomainChange>& domchgstack,
                            const std::vector<HighsInt>& branchingPositions);

  bool propagate(bool dualflag = false);

  double getColLowerPos(HighsInt col, HighsInt stackpos, HighsInt& pos) const;

  double getColUpperPos(HighsInt col, HighsInt stackpos, HighsInt& pos) const;

  void conflictAnalysis(HighsConflictPool& conflictPool);

  void conflictAnalysis(const HighsInt* proofinds, const double* proofvals,
                        HighsInt prooflen, double proofrhs,
                        HighsConflictPool& conflictPool);

  void conflictAnalyzeReconvergence(const HighsDomainChange& domchg,
                                    const HighsInt* proofinds,
                                    const double* proofvals, HighsInt prooflen,
                                    double proofrhs,
                                    HighsConflictPool& conflictPool);

  void tightenCoefficients(HighsInt* inds, double* vals, HighsInt len,
                           double& rhs) const;

  double getMinActivity(HighsInt row) const {
    return activitymininf_[row] == 0 ? double(activitymin_[row]) : -kHighsInf;
  }

  double getMaxActivity(HighsInt row) const {
    return activitymaxinf_[row] == 0 ? double(activitymax_[row]) : kHighsInf;
  }

  double getMinCutActivity(const HighsCutPool& cutpool, HighsInt cut) const;

  bool isBinary(HighsInt col) const {
    return mipsolver->isColIntegral(col) && col_lower_[col] == 0.0 &&
           col_upper_[col] == 1.0;
  }

  bool isGlobalBinary(HighsInt col) const {
    return mipsolver->isColIntegral(col) &&
           mipsolver->model_->col_lower_[col] == 0.0 &&
           mipsolver->model_->col_upper_[col] == 1.0;
  }

  HighsVarType variableType(HighsInt col) const {
    return mipsolver->variableType(col);
  }

  bool isFixed(HighsInt col) const {
    return col_lower_[col] == col_upper_[col];
  }

  bool isFixing(const HighsDomainChange& domchg) const;

  HighsDomainChange flip(const HighsDomainChange& domchg) const;

  double feastol() const;

  HighsInt numModelNonzeros() const { return mipsolver->numNonzero(); }

  bool inSubmip() const { return mipsolver->submip; }

  void clearRedundantRows() { redundantRows_.clear(); };

  const HighsHashTable<HighsInt>& getRedundantRows() const {
    return redundantRows_;
  };

  double getRedundantRowValue(HighsInt row) const;

  void setRecordRedundantRows(bool val) { recordRedundantRows_ = val; };

  bool isRedundantRow(HighsInt row) const;
};

#endif
