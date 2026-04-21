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
    components.save_feedback({"path_a": 1, "path_b": 0}, study / "path_feedback.parquet", "alice", "my_study")
    components.save_feedback({"rxn_a": 1}, study / "reaction_feedback.parquet", "alice", "my_study")

    path_fb, rxn_fb = components.load_user_feedback("alice", "my_study")
    assert path_fb == {"path_a": 1, "path_b": 0}
    assert rxn_fb == {"rxn_a": 1}


def test_upsert_updates_existing_row(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", "alice", "my_study")
    components.save_feedback({"path_a": 0}, study / "path_feedback.parquet", "alice", "my_study")

    path_fb, _ = components.load_user_feedback("alice", "my_study")
    assert path_fb == {"path_a": 0}


def test_usernames_aggregated_across_tables(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", "alice", "my_study")
    components.save_feedback({"rxn_a": 1}, study / "reaction_feedback.parquet", "bob", "my_study")

    assert components.get_existing_usernames() == ["alice", "bob"]


def test_guest_is_noop(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", "guest", "my_study")
    assert components.load_user_feedback("guest", "my_study") == ({}, {})
    assert components.get_existing_usernames() == []


def test_feedback_isolated_per_user(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", "alice", "my_study")
    components.save_feedback({"path_a": 0}, study / "path_feedback.parquet", "bob", "my_study")

    alice_fb, _ = components.load_user_feedback("alice", "my_study")
    bob_fb, _ = components.load_user_feedback("bob", "my_study")
    assert alice_fb == {"path_a": 1}
    assert bob_fb == {"path_a": 0}


def test_feedback_isolated_per_study(feedback_env):
    _, study = feedback_env
    components.save_feedback({"path_a": 1}, study / "path_feedback.parquet", "alice", "study_a")
    components.save_feedback({"path_a": 0}, study / "path_feedback.parquet", "alice", "study_b")

    a_fb, _ = components.load_user_feedback("alice", "study_a")
    b_fb, _ = components.load_user_feedback("alice", "study_b")
    assert a_fb == {"path_a": 1}
    assert b_fb == {"path_a": 0}

    c_fb, _ = components.load_user_feedback("alice", "other_study")
    assert c_fb == {}


def test_unknown_filepath_raises(feedback_env):
    _, study = feedback_env
    with pytest.raises(ValueError):
        components.save_feedback({"x": 1}, study / "bogus.parquet", "alice", "my_study")


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

    path_fb, _ = components.load_user_feedback("alice", "my_study")
    assert path_fb == {"path_a": 1, "path_b": 0}

    other, _ = components.load_user_feedback("alice", "some_other_study")
    assert other == {}


def test_migration_skipped_when_db_nonempty(feedback_env):
    _, study_root = feedback_env
    components.save_feedback({"path_a": 1}, study_root / "path_feedback.parquet", "alice", "my_study")

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

    path_fb, _ = components.load_user_feedback("alice", "my_study")
    assert path_fb == {"path_a": 1}
