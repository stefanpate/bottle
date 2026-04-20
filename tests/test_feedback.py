import polars as pl
import pytest

from path_viewer import components


@pytest.fixture
def feedback_env(tmp_path, monkeypatch):
    fb_root = tmp_path / "feedback"
    study_root = tmp_path / "processed"
    fb_root.mkdir()
    study_root.mkdir()
    monkeypatch.setenv("FEEDBACK_ROOT", str(fb_root))
    monkeypatch.setenv("CASP_STUDY_ROOT", str(study_root))
    return fb_root, study_root


def test_save_and_load_roundtrip(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1, "path_b": 0}, study / "path_feedback.parquet", None, "alice")
    components.save_feedback({"rxn_a": 1}, study / "reaction_feedback.parquet", None, "alice")

    path_fb, rxn_fb = components.load_user_feedback("alice", study)
    assert path_fb == {"path_a": 1, "path_b": 0}
    assert rxn_fb == {"rxn_a": 1}


def test_upsert_updates_existing_row(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", None, "alice")
    components.save_feedback({"path_a": 0}, study / "path_feedback.parquet", None, "alice")

    path_fb, _ = components.load_user_feedback("alice", study)
    assert path_fb == {"path_a": 0}


def test_usernames_aggregated_across_tables(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", None, "alice")
    components.save_feedback({"rxn_a": 1}, study / "reaction_feedback.parquet", None, "bob")

    assert components.get_existing_usernames() == ["alice", "bob"]


def test_guest_is_noop(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", None, "guest")
    assert components.load_user_feedback("guest", study) == ({}, {})
    assert components.get_existing_usernames() == []


def test_feedback_isolated_per_user(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", None, "alice")
    components.save_feedback({"path_a": 0}, study / "path_feedback.parquet", None, "bob")

    alice_fb, _ = components.load_user_feedback("alice", study)
    bob_fb, _ = components.load_user_feedback("bob", study)
    assert alice_fb == {"path_a": 1}
    assert bob_fb == {"path_a": 0}


def test_unknown_filepath_raises(feedback_env):
    _, study = feedback_env
    with pytest.raises(ValueError):
        components.save_feedback({"x": 1}, study / "bogus.parquet", None, "alice")


def test_legacy_parquet_migration(feedback_env):
    _, study_root = feedback_env
    study = study_root / "my_study"
    study.mkdir()
    legacy = pl.DataFrame(
        {
            "username": ["alice", "alice"],
            "id": ["path_a", "path_b"],
            "feedback": [1, 0],
            "date": ["2026-04-20", "2026-04-20"],
            "time": ["12:00:00", "12:01:00"],
        }
    )
    legacy.write_parquet(study / "path_feedback.parquet")

    path_fb, _ = components.load_user_feedback("alice", study)
    assert path_fb == {"path_a": 1, "path_b": 0}


def test_migration_skipped_when_db_nonempty(feedback_env):
    _, study_root = feedback_env
    components.save_feedback({"path_a": 1}, study_root / "path_feedback.parquet", None, "alice")

    study = study_root / "my_study"
    study.mkdir()
    legacy = pl.DataFrame(
        {
            "username": ["alice"],
            "id": ["path_a"],
            "feedback": [0],
            "date": ["2026-04-20"],
            "time": ["12:00:00"],
        }
    )
    legacy.write_parquet(study / "path_feedback.parquet")

    path_fb, _ = components.load_user_feedback("alice", study)
    assert path_fb == {"path_a": 1}
